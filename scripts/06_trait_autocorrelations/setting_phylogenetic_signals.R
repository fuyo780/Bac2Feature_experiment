library(ape)
library(castor)
library(dplyr)
library(ggplot2)
library(reshape)


# import tree and trait data
tree.file <- "../../data/ref_bac2feature/phylogeny/phylogeny.tre"
trait.file <- "../../data/ref_bac2feature/trait_bac2feature.tsv"

tree <- ape::read.tree(tree.file)

traits <- read.table(trait.file, sep="\t", comment.char="", quote = "", header=T)
traits$species_tax_id <- as.character(traits$species_tax_id)
rownames(traits) <- traits$species_tax_id

# Auto-correlation result
autocorrelation <- read.table("../../data/trait_autocorrelations/autocorrelation.tsv", sep="\t", comment.char="", quote = "", header=F, col.names=c("trait", "lower", "mean", "upper", "phylogenetic_distance"))

continuous_traits <- c("cell_diameter", "cell_length", "doubling_h", "genome_size", "gc_content", "coding_genes", "optimum_tmp", "optimum_ph", "growth_tmp", "rRNA16S_genes", "tRNA_genes")

categorical_traits <- c('gram_stain','sporulation','motility','range_salinity','facultative_respiration','anaerobic_respiration','aerobic_respiration','mesophilic_range_tmp','thermophilic_range_tmp','psychrophilic_range_tmp','bacillus_cell_shape','coccus_cell_shape','filament_cell_shape','coccobacillus_cell_shape','vibrio_cell_shape','spiral_cell_shape')

min_phylodistance <- 0.0001
max_phylodistance <- 2

all_traits <- c(continuous_traits, categorical_traits)

# Calculate the threshold values of phylogenetic distances.
ac_copy <- autocorrelation
# Fill NA in confidence interval
fill_na <- function(df, col1, col2) { # Fill col1 with col2 if col1 is NA
  return(ifelse(is.na(df[col1]), df[col2], df[col1]))
}
ac_copy[,"lower"] <- as.numeric(apply(ac_copy, 1, fill_na, "lower", "mean"))
ac_copy[,"upper"] <- as.numeric(apply(ac_copy, 1, fill_na, "upper", "mean"))

threshold_phylodistance_strict <- c()
threshold_correlation <- 0.5
# Take one of the previous values of the phylogenetic distance when the correlation falls below 0 for the first time
for (t in all_traits) {
  ac_trait <- ac_copy %>% filter(trait == t)
  under_threshold_idx <- which(ac_trait["lower"] <= threshold_correlation)[1]
  if (min(ac_trait$lower) > threshold_correlation) {
    threshold_phylodistance_strict[t] <- max_phylodistance
  } else if (max(ac_trait$lower) < threshold_correlation) {
    threshold_phylodistance_strict[t] <- 0
  } else if (under_threshold_idx == 1) {
    threshold_phylodistance_strict[t] <- 0
  } else {
    threshold_phylodistance_strict[t] <- ac_trait[under_threshold_idx-1, "phylogenetic_distance"]
  }
}

threshold_phylodistance_loose <- c()
threshold_correlation <- 0.0
# Take one of the previous values of the phylogenetic distance when the correlation falls below 0 for the first time
for (t in all_traits) {
  ac_trait <- ac_copy %>% filter(trait == t)
  under_threshold_idx <- which(ac_trait["lower"] <= threshold_correlation)[1]
  if (min(ac_trait$lower) > threshold_correlation) {
    threshold_phylodistance_loose[t] <- max_phylodistance
  } else if (max(ac_trait$lower) < threshold_correlation) {
    threshold_phylodistance_loose[t] <- 0
  } else if (under_threshold_idx == 1) {
    threshold_phylodistance_loose[t] <- 0
  } else {
    threshold_phylodistance_loose[t] <- ac_trait[under_threshold_idx-1, "phylogenetic_distance"]
  }
}

n <- length(all_traits)
res_df <- data.frame(
  trait=all_traits,
  cor_0.5=threshold_phylodistance_strict,
  cor_0=threshold_phylodistance_loose
)

write.table(res_df, "../../data/trait_autocorrelations/threshold_phylodistance.tsv", sep = "\t", row.names = F, col.names = T, quote = F)
