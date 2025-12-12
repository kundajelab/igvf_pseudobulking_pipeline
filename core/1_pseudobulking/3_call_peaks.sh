#!/bin/bash

callpeak () {
	dataset=${1}
	datadir=${2}
	chr_order=${3}
	blacklist=${4}

	p1_dir="${datadir}/pseudobulked_pseudorep1"
	p2_dir="${datadir}/pseudobulked_pseudorep2"
	pT_dir="${datadir}/pseudobulked_pseudorepT"
	out_dir="${datadir}/peaks"

	: '
	if [ -f "${out_dir}/${dataset}_peaks_overlap_filtered.narrowPeak" ]; then
		echo "${dataset} already completed! skipping..."
		return 1
	fi
	'

	echo -e "\t\t- ${dataset} calling pseudorep1 peaks..."
	macs3 callpeak -t ${p1_dir}/${dataset}.tsv -f BED -n ${dataset}-pseudoreplicate1 -g hs --outdir ${out_dir} -p 0.01 --shift -75 --extsize 150 --nomodel -B --SPMR --keep-dup all --call-summits &> ${out_dir}/log-${dataset}-pseudorep1.txt
	echo -e "\t\t- ${dataset} calling pseudorep2 peaks..."
	macs3 callpeak -t ${p2_dir}/${dataset}.tsv -f BED -n ${dataset}-pseudoreplicate2 -g hs --outdir ${out_dir} -p 0.01 --shift -75 --extsize 150 --nomodel -B --SPMR --keep-dup all --call-summits &> ${out_dir}/log-${dataset}-pseudorep2.txt
	echo -e "\t\t- ${dataset} calling pseudorepT peaks..."
	macs3 callpeak -t ${pT_dir}/${dataset}.tsv -f BED -n ${dataset}-pseudoreplicateT -g hs --outdir ${out_dir} -p 0.01 --shift -75 --extsize 150 --nomodel -B --SPMR --keep-dup all --call-summits &> ${out_dir}/log-${dataset}-pseudorepT.txt

	p1_in=${out_dir}/${dataset}-pseudoreplicate1_peaks.narrowPeak
	p2_in=${out_dir}/${dataset}-pseudoreplicate2_peaks.narrowPeak
	pT_in=${out_dir}/${dataset}-pseudoreplicateT_peaks.narrowPeak
	p1_out=${out_dir}/${dataset}-pseudoreplicate1_peaks_top.narrowPeak
	p2_out=${out_dir}/${dataset}-pseudoreplicate2_peaks_top.narrowPeak
	pT_out=${out_dir}/${dataset}-pseudoreplicateT_peaks_top.narrowPeak
	npeaks=300000

	echo -e "\t\t- ${dataset} getting top peaks..."
	sort -k 8gr,8gr ${p1_in} | head -n ${npeaks} | bedtools sort -i stdin -g ${chr_order} > ${p1_out}
	sort -k 8gr,8gr ${p2_in} | head -n ${npeaks} | bedtools sort -i stdin -g ${chr_order} > ${p2_out}
	sort -k 8gr,8gr ${pT_in} | head -n ${npeaks} | bedtools sort -i stdin -g ${chr_order} > ${pT_out}

	echo -e "\t\t- ${dataset} intersecting peaks"
	min_overlap=0.5
	overlap_output=${out_dir}/${dataset}-peaks_overlap.narrowPeak
	bedtools intersect -u -a ${pT_out} -b ${p1_out} -g ${chr_order} -f ${min_overlap} -F ${min_overlap} -e -sorted | bedtools intersect -u -a stdin -b ${p2_out} -g ${chr_order} -f ${min_overlap} -F ${min_overlap} -e -sorted > ${overlap_output}

	echo -e "\t\t- ${dataset} filtering blacklist peaks"
	filtered_output=${out_dir}/${dataset}-peaks_overlap_filtered.narrowPeak
	bedtools intersect -v -a ${overlap_output} -b ${blacklist} > ${filtered_output}

	echo -e "\t\t- ${dataset} making p-value bedgraphs"
	macs3 bdgcmp -m ppois -t ${out_dir}/${dataset}-pseudoreplicate1_treat_pileup.bdg -c ${out_dir}/${dataset}-pseudoreplicate1_control_lambda.bdg -o ${out_dir}/${dataset}-pseudoreplicate1_ppois.bdg
	macs3 bdgcmp -m ppois -t ${out_dir}/${dataset}-pseudoreplicate2_treat_pileup.bdg -c ${out_dir}/${dataset}-pseudoreplicate2_control_lambda.bdg -o ${out_dir}/${dataset}-pseudoreplicate2_ppois.bdg
	macs3 bdgcmp -m ppois -t ${out_dir}/${dataset}-pseudoreplicateT_treat_pileup.bdg -c ${out_dir}/${dataset}-pseudoreplicateT_control_lambda.bdg -o ${out_dir}/${dataset}-pseudoreplicateT_ppois.bdg
	bedClip ${out_dir}/${dataset}-pseudoreplicate1_ppois.bdg ${chr_order} ${out_dir}/${dataset}-pseudoreplicate1_ppois_clipped.bdg
	bedClip ${out_dir}/${dataset}-pseudoreplicate2_ppois.bdg ${chr_order} ${out_dir}/${dataset}-pseudoreplicate2_ppois_clipped.bdg
	bedClip ${out_dir}/${dataset}-pseudoreplicateT_ppois.bdg ${chr_order} ${out_dir}/${dataset}-pseudoreplicateT_ppois_clipped.bdg
	echo -e "\t\t- ${dataset} combining p-value bedgraphs"
	macs3 cmbreps -m fisher -i ${out_dir}/${dataset}-pseudoreplicate1_ppois_clipped.bdg ${out_dir}/${dataset}-pseudoreplicate2_ppois_clipped.bdg ${out_dir}/${dataset}-pseudoreplicateT_ppois_clipped.bdg -o ${out_dir}/${dataset}-combined_ppois.bdg
	tail -n +2 ${out_dir}/${dataset}-combined_ppois.bdg | bedtools sort -i stdin -g ${chr_order} > ${out_dir}/${dataset}-combined_ppois_sorted.bdg
	echo -e "\t\t- ${dataset} converting p-value bedgraphs to bigwigs"
	bedGraphToBigWig ${out_dir}/${dataset}-combined_ppois_sorted.bdg ${chr_order} ${out_dir}/${dataset}-pval.bw

	echo -e "\t\t- ${dataset} making raw insertion bigwig"
	bedtools genomecov -i ${pT_dir}/${dataset}.tsv -g ${chr_order} -bg > ${out_dir}/${dataset}-raw_insertions.bdg
	bedGraphToBigWig ${out_dir}/${dataset}-raw_insertions.bdg ${chr_order} ${out_dir}/${dataset}-raw_insertions.bw

	echo -e "\t\t- ${dataset} computing per cell FRiP"
	frags_dir="${datadir}/pseudobulked_fragments"
	cut -f4 "${frags_dir}/${dataset}-sorted.tsv" | sort | uniq -c | awk '{print $2, $1}' > "${out_dir}/${dataset}-fragments_per_cell.txt"
	bedtools intersect -a  "${frags_dir}/${dataset}-sorted.tsv" -b ${filtered_output} -u > "${out_dir}/${dataset}-fragments_in_peaks.tsv"
	cut -f4 "${out_dir}/${dataset}-fragments_in_peaks.tsv" | sort | uniq -c | awk '{print $2, $1}' > "${out_dir}/${dataset}-fragments_in_peaks_per_cell.txt"
	join -a 1 -1 1 -2 1 <(sort "${out_dir}/${dataset}-fragments_per_cell.txt") <(sort "${out_dir}/${dataset}-fragments_in_peaks_per_cell.txt") | awk '{if($3=="") $3=0; frip=$3/$2; print $1, frip}' > "${out_dir}/${dataset}-frip_per_cell.txt"

	echo -e "\t\t${dataset} completed!"

}
export -f callpeak

datadir="${1}"
chr_order_file="${2}"
blacklist_file="${3}"
parallel="${4}"

frags_dir="${datadir}/pseudobulked_fragments"

ls ${frags_dir} | grep sorted | sed 's/-[^-]*$//' | sort | uniq | parallel --linebuffer -j ${parallel} callpeak {} ${datadir} ${chr_order_file} ${blacklist_file}
