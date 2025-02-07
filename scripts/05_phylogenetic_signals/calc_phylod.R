library(castor)
library(dplyr)
library(ggplot2)
library(picante)
library(caper)

# Phylogenetic tree
tree_path <- "../../data/ref_bac2feature/phylogeny/phylogeny.tre"
tree <- read.tree(tree_path)
# Trait
trait_path <- "../../data/ref_bac2feature/trait_bac2feature.tsv"
trait <- read.table(trait_path, sep="\t", comment.char="", quote = "", header=T)
trait$species_tax_id <- as.character(trait$species_tax_id)
rownames(trait) <- trait$species_tax_id

# Use modified phylod function
source("./modified_phylod.R")

traits_binned <- trait[c(c(1), c(13:28))]
ct <- colnames(traits_binned)[-1]

# Calculation of D metrics based on original tree
l_phylod <- c()
for (c in ct) {
  known_trait_inds <- which(!is.na(traits_binned[,c]))
  known_trait_cols <- c("species_tax_id", c)
  known_traits <- traits_binned[known_trait_inds, known_trait_cols]
  subtree <- get_subtree_with_tips(tree=tree,
                                   only_tips=known_traits$species_tax_id
                                  )$subtree
  phylod <- phylo.d2(data=known_traits,
                     phy=subtree,
                     names.col=species_tax_id,
                     binvar=c,
                     permut=1000,
                     rnd.bias=NULL
                    )
  l_phylod <- c(l_phylod, c(phylod$DEstimate))
}

# Save the results
names(l_phylod) <- ct
df_phylod <- data.frame(l_phylod)
write.table(df_phylod,
            file="../../data/phylogenetic_signals/phylod.tsv",
            sep = "\t",
            na="",
            col.names=F,
            quote=F
            )
