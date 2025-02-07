# This script has to be called in b2f environment

library(ape)
library(optparse)

option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL,
                help="Input tree in newick format", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="Output subtree in newick", metavar="character"),
    make_option(c("-k", "--keep"), type="list", default=NULL,
                help="File including character vector listing tip names to keep (Note even interger is classified as tip names)", metavar="character"),
    make_option(c("-p", "--prune"), type="list", default=NULL,
                help="File including character vector listing tip names to prune (Note even interger is classified as tip names)", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$output)) {
    print_help(opt_parser)
    stop("At least one argument must be supplied (input and output file)\n", call.=FALSE)
}
if (is.null(opt$keep) && is.null(opt$prune)) {
    print_help(opt_parser)
    stop("At least one argument must be supplied (keep or prune argument)\n", call.=FALSE)
}
if (!is.null(opt$keep) && !is.null(opt$prune)) {
    print_help(opt_parser)
    stop("Only one argument must be supplied (keep or prune argument)\n", call.=FALSE)
}

tree = read.tree(opt$input)

keep_tips = NULL
prune_tips = NULL
pruned_tree = NULL

if (!is.null(opt$keep)) {
    df = read.table(opt$keep)
    keep_tips = as.character(df[,1]) # input is a vector of tip name
    pruned_tree = keep.tip(tree, tip=keep_tips)
}
if (!is.null(opt$prune)) {
    df = read.table(opt$prune)
    prune_tips = as.character(df[,1]) # input is a vector of tip name
    pruned_tree = drop.tip(tree, tip=prune_tips)
}

write.tree(phy=pruned_tree,
           file=opt$output,
           append=FALSE,
           digits=10,
           tree.names=FALSE)
