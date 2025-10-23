#!/bin/bash

catsort () {
	dataset=${1}
	datadir=${2}

	in_base="${datadir}/separated"
	out_base="${datadir}/pseudobulked"

	echo -e "\t- ${dataset} concatenating across samples..."
	find "${in_base}_fragments" -name "${dataset}*.tsv" -exec cat {} + > "${out_base}_fragments/${dataset}.tsv"
	find "${in_base}_pseudorepT" -name "${dataset}*.tsv" -exec cat {} + > "${out_base}_pseudorepT/${dataset}.tsv"
	find "${in_base}_pseudorep1" -name "${dataset}*.tsv" -exec cat {} + > "${out_base}_pseudorep1/${dataset}.tsv"
	find "${in_base}_pseudorep2" -name "${dataset}*.tsv" -exec cat {} + > "${out_base}_pseudorep2/${dataset}.tsv"

	echo -e "\t- ${dataset} sorting fragments..."
	sort -k 1,1 -k 2,2n -S 10% --parallel=4 ${out_base}_fragments/${dataset}.tsv -o ${out_base}_fragments/${dataset}-sorted.tsv
	echo -e "\t- ${dataset} sorting pseudorepT..."
	sort -k 1,1 -k 2,2n -S 10% --parallel=4 ${out_base}_pseudorepT/${dataset}.tsv -o ${out_base}_pseudorepT/${dataset}-sorted.tsv
	echo -e "\t- ${dataset} sorting pseudorep1..."
	sort -k 1,1 -k 2,2n -S 10% --parallel=4 ${out_base}_pseudorep1/${dataset}.tsv -o ${out_base}_pseudorep1/${dataset}-sorted.tsv
	echo -e "\t- ${dataset} sorting pseudorep2..."
	sort -k 1,1 -k 2,2n -S 10% --parallel=4 ${out_base}_pseudorep2/${dataset}.tsv -o ${out_base}_pseudorep2/${dataset}-sorted.tsv

	echo -e "\t- done ${dataset}"
}
export -f catsort

datadir="${1}"
parallel="${2}"

frags_dir="${datadir}/separated_fragments"

ls ${frags_dir} | cut -d "-" -f 1 | sort | uniq | parallel --linebuffer -j ${parallel} catsort {} ${datadir}

