library(ape)
library(castor)
library(dplyr)
library(ggplot2)
library(reshape)

library(phylobase)
library(phylosignal)

# Read in command-line arguments.
Args <- commandArgs(TRUE)

# import tree and trait data
tree.file <- "../../data/ref_bac2feature/phylogeny/phylogeny.tre"
trait.file <- "../../data/ref_bac2feature/trait_bac2feature.tsv"

tree <- ape::read.tree(tree.file)

traits <- read.table(trait.file, sep="\t", comment.char="", quote = "", header=T)
traits$species_tax_id <- as.character(traits$species_tax_id)
rownames(traits) <- traits$species_tax_id

# func
get_distance_df <- function(tree, traits, t) {
  # Prune tree
  known_spp <- traits[!is.na(traits[,t]), "species_tax_id"]
  subtree <- castor::get_subtree_with_tips(tree = tree, only_tips = known_spp)$subtree
  subtrait <- traits[known_spp, t]
  Ntips <- length(subtree$tip.label)

  # Calculate distance matrix
  dist_all <- castor::get_all_pairwise_distances(tree = subtree)
  dist_tip <- dist_all[1:Ntips, 1:Ntips]
  rownames(dist_tip) <- subtree$tip.label
  colnames(dist_tip) <- subtree$tip.label
  dist_df <- melt.matrix(data = dist_tip, as.is = TRUE)
  colnames(dist_df) <- c("tip1", "tip2", "dist")
  dist_df$tip1 <- as.character(dist_df$tip1)
  dist_df$tip2 <- as.character(dist_df$tip2)
  return(dist_df)
}

calc_trait_corr <- function(all_pairs_in_bin, t, lw, up, Npairs) {
  # Subsample from all species pairs
  sub_pairs_in_bin <- sample_n(all_pairs_in_bin, min(Npairs, dim(all_pairs_in_bin)[1]), replace = F)
  # Merge trait info
  xtraits <- traits
  colnames(xtraits) <- paste0("x.", colnames(xtraits))
  sub_pairs_x <- left_join(x = sub_pairs_in_bin, y = xtraits, by = c("tip1" = "x.species_tax_id"))
  ytraits <- traits
  colnames(ytraits) <- paste0("y.", colnames(ytraits))
  sub_pairs_xy <- left_join(x = sub_pairs_x, y = ytraits, by = c("tip2" = "y.species_tax_id"))
  # Calculate correlation
  x <- sub_pairs_xy[, paste0("x.", t)]
  y <- sub_pairs_xy[, paste0("y.", t)]
  return(cor(x, y, method = "pearson"))
}

set.seed(123)

min_phylodistance <- 0.0001
max_phylodistance <- 2
bin_width <- 0.05
Npairs <- 1000
Niter <- 100
outlier_per <- 0.05

trait_cols <- colnames(traits)[-1]

bins <- seq(0, max_phylodistance, bin_width)
bins[1] <- min_phylodistance

for (t in trait_cols) {
  dist_df <- get_distance_df(tree = tree, traits = traits, t = t)

  cor_quantile_matrix <- c()
  for (i in seq(1, length(bins)-1)) {
    lw <- bins[i]
    up <- bins[i+1]
    cor_list <- c() # correlation list
    rm_ci_flag <- FALSE
    all_pairs_in_bin <- dist_df[dist_df$dist > lw & dist_df$dist < up, c("tip1", "tip2")] # distance matrix
    if (dim(all_pairs_in_bin)[1] < 1) { # No species in bin: skip
      cor_list <- c(NaN)
    } else if (dim(all_pairs_in_bin)[1] <= Npairs) { # Npairs or less species in bin: calculate correlation once
      cor_list <- c(calc_trait_corr(all_pairs_in_bin = all_pairs_in_bin, t = t, lw = lw, up = up, Npairs = Npairs))
      rm_ci_flag <- TRUE
    } else {
      for (i in seq(Niter)) { # more than Npairs species in bin: calculate correlation Niter times
        cor <- calc_trait_corr(all_pairs_in_bin = all_pairs_in_bin, t = t, lw = lw, up = up, Npairs = Npairs)
        cor_list <- c(cor_list, c(cor))
      }
    }
    quantile_list <- quantile(cor_list, c(outlier_per / 2, 0.50, 1 - outlier_per / 2), na.rm = T)
    if (rm_ci_flag) {
      quantile_list[[1]] <- NaN
      quantile_list[[3]] <- NaN
    }
    quantile_list["phylogenetic_distance"] <- up
    cor_quantile_matrix <- rbind(cor_quantile_matrix, quantile_list)
  }

  rownames(cor_quantile_matrix) <- NULL
  colnames(cor_quantile_matrix) <- c("lw", "mean", "up", "phylogenetic_distance")

  cor_quantile_df <- as.data.frame(cor_quantile_matrix)
  cor_quantile_df["trait"] <- t
  cor_quantile_df <- cor_quantile_df[,c(5,1,2,3,4)]

  write.table(cor_quantile_df, "../../data/trait_autocorrelations/autocorrelation.tsv", sep = "\t", row.names = F, col.names = F, quote = F, append = T)
}
