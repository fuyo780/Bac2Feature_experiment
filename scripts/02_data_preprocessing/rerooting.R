library(ape)
library(castor)

b2f_tree <- read.tree("../../data/intermediate_dir/phylogeny.raxml.bestTree")
rerooted_b2f_tree <- castor::root_at_midpoint(b2f_tree, castor::find_root(b2f_tree))

write_tree(rerooted_b2f_tree, "../../data/intermediate_dir/phylogeny.raxml.bestTree.rerooted")
