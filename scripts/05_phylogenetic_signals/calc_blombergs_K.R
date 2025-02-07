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

# It takes a lot of time: ~120min
blombergs_K <- multiPhylosignal(trait[tree$tip.label, c(2:12)], tree)
write.table(blombergs_K,
            file="../../data/phylogenetic_signals/blombergs_K.tsv",
            sep = "\t",
            na="",
            row.names=F,
            quote=F
            )
