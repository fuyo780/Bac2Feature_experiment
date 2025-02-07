#!/bin/bash

while getopts t:a:n:o: option
do
  case $option in
    t ) tree_path=$OPTARG;;
    a ) msa_path=$OPTARG;;
    n ) node_id_path=$OPTARG;;
    o ) outdir_path=$OPTARG;;
    \?) echo "This is unexpected option." 1>&2
        exit 1
  esac
done

input_file_path=($tree_path $msa_path $node_id_path)
for f in ${input_file_path[@]}
do
  if [ ! -f $f ];then
    echo "Input file doesn't exits: ${f}"
    exit 1
  fi
done

lib_dir=$(dirname $0)
out_prefix=$(basename $outdir_path)
out_prefix_path="${outdir_path%/}/${out_prefix}"
tmp_start_tree_path="${out_prefix_path}.startTree"
out_tree_path="${out_prefix_path}.tre"
out_msa_path="${out_prefix_path}.fasta"
out_hmm_path="${out_prefix_path}.hmm"
out_model_path="${out_prefix_path}.model"

if [ ! -d $outdir_path ];then
  mkdir $outdir_path
fi

# Prune tree according to node_id list
Rscript $lib_dir/prune_tree.R -i $tree_path -o $tmp_start_tree_path -k $node_id_path

# Extract MSA according to node_id list
seqkit grep -nf $node_id_path $msa_path 1> $out_msa_path 2> /dev/null

# Predict HMM profile
hmmbuild $out_hmm_path $out_msa_path > /dev/null 2>&1

# Optimize by Raxml-ng
raxml-ng --evaluate --msa $out_msa_path --tree $tmp_start_tree_path --prefix $out_prefix_path --model GTR+G+F > /dev/null 2>&1
# Use the best tree and model
cp "${out_prefix_path}.raxml.bestTree" $out_tree_path
cp "${out_prefix_path}.raxml.bestModel" $out_model_path
rm ${out_prefix_path}.raxml.* $tmp_start_tree_path
