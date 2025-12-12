#!/bin/bash

catsort () {
	dataset=${1}
	datadir=${2}
	chrom_sizes=${3}

	in_base="${datadir}/separated"
	out_base="${datadir}/pseudobulked"

	echo -e "\t\t- ${dataset} concatenating across samples..."
	find "${in_base}_fragments" -name "${dataset}*.tsv" -exec cat {} + > "${out_base}_fragments/${dataset}.tsv"
	find "${in_base}_pseudorepT" -name "${dataset}*.tsv" -exec cat {} + > "${out_base}_pseudorepT/${dataset}.tsv"
	find "${in_base}_pseudorep1" -name "${dataset}*.tsv" -exec cat {} + > "${out_base}_pseudorep1/${dataset}.tsv"
	find "${in_base}_pseudorep2" -name "${dataset}*.tsv" -exec cat {} + > "${out_base}_pseudorep2/${dataset}.tsv"

	echo -e "\t\t- ${dataset} sorting fragments..."
	bedtools sort -i ${out_base}_fragments/${dataset}.tsv -g ${chrom_sizes} > ${out_base}_fragments/${dataset}-sorted.tsv
	# echo -e "\t\t- ${dataset} sorting pseudorepT..."
	# bedtools sort -i ${out_base}_pseudorepT/${dataset}.tsv -g ${chrom_sizes} > ${out_base}_pseudorepT/${dataset}-sorted.tsv
	# echo -e "\t\t- ${dataset} sorting pseudorep1..."
	# bedtools sort -i ${out_base}_pseudorep1/${dataset}.tsv -g ${chrom_sizes} > ${out_base}_pseudorep1/${dataset}-sorted.tsv
	# echo -e "\t\t- ${dataset} sorting pseudorep2..."
	# bedtools sort -i ${out_base}_pseudorep2/${dataset}.tsv -g ${chrom_sizes} > ${out_base}_pseudorep2/${dataset}-sorted.tsv

	echo -e "\t\t- done ${dataset}"
}
export -f catsort

datadir="${1}"
parallel="${2}"
chrom_sizes="${3}"

frags_dir="${datadir}/separated_fragments"

ls ${frags_dir} | sed 's/-[^-]*$//' | sort | uniq | parallel --linebuffer -j ${parallel} catsort {} ${datadir} ${chrom_sizes}

