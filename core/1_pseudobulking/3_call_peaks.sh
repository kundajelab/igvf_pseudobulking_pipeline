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

	echo -e "\t- ${dataset} calling pseudorep1 peaks..." &
	macs3 callpeak -t ${p1_dir}/${dataset}-sorted.tsv -f BED -n ${dataset}-pseudoreplicate1 -g hs --outdir ${out_dir} -p 0.01 --shift -75 --extsize 150 --nomodel -B --SPMR --keep-dup all --call-summits &> ${out_dir}/log-${dataset}-pseudorep1.txt &
	echo -e "\t- ${dataset} calling pseudorep2 peaks..." &
	macs3 callpeak -t ${p2_dir}/${dataset}-sorted.tsv -f BED -n ${dataset}-pseudoreplicate2 -g hs --outdir ${out_dir} -p 0.01 --shift -75 --extsize 150 --nomodel -B --SPMR --keep-dup all --call-summits &> ${out_dir}/log-${dataset}-pseudorep2.txt &
	echo -e "\t- ${dataset} calling pseudorepT peaks..." &
	macs3 callpeak -t ${pT_dir}/${dataset}-sorted.tsv -f BED -n ${dataset}-pseudoreplicateT -g hs --outdir ${out_dir} -p 0.01 --shift -75 --extsize 150 --nomodel -B --SPMR --keep-dup all --call-summits &> ${out_dir}/log-${dataset}-pseudorepT.txt &
	wait
	echo -e "\t- ${dataset} finished waiting"

	p1_in=${out_dir}/${dataset}-pseudoreplicate1_peaks.narrowPeak
	p2_in=${out_dir}/${dataset}-pseudoreplicate2_peaks.narrowPeak
	pT_in=${out_dir}/${dataset}-pseudoreplicateT_peaks.narrowPeak
	p1_out=${out_dir}/${dataset}-pseudoreplicate1_peaks_top.narrowPeak
	p2_out=${out_dir}/${dataset}-pseudoreplicate2_peaks_top.narrowPeak
	pT_out=${out_dir}/${dataset}-pseudoreplicateT_peaks_top.narrowPeak
	npeaks=300000

	echo -e "\t- ${dataset} getting top peaks..."
        sort -k 8gr,8gr ${p1_in} | head -n ${npeaks} | sort -k 1,1 -k2,2n > ${p1_out}
        sort -k 8gr,8gr ${p2_in} | head -n ${npeaks} | sort -k 1,1 -k2,2n > ${p2_out}
        sort -k 8gr,8gr ${pT_in} | head -n ${npeaks} | sort -k 1,1 -k2,2n > ${pT_out}

	echo -e "\t- ${dataset} intersecting peaks"
        min_overlap=0.5
        overlap_output=${out_dir}/${dataset}-peaks_overlap.narrowPeak
        bedtools intersect -u -a ${pT_out} -b ${p1_out} -g ${chr_order} -f ${min_overlap} -F ${min_overlap} -e -sorted | bedtools intersect -u -a stdin -b ${p2_out} -g ${chr_order} -f ${min_overlap} -F ${min_overlap} -e -sorted > ${overlap_output}

	echo -e "\t- ${dataset} filtering blacklist peaks"
        filtered_output=${out_dir}/${dataset}-peaks_overlap_filtered.narrowPeak
        bedtools intersect -v -a ${overlap_output} -b ${blacklist} > ${filtered_output}

	echo -e "\t- ${dataset} making p-value bedgraphs"
	macs3 bdgcmp -m ppois -t ${out_dir}/${dataset}-pseudoreplicate1_treat_pileup.bdg -c ${out_dir}/${dataset}-pseudoreplicate1_control_lambda.bdg -o ${out_dir}/${dataset}-pseudoreplicate1_ppois.bdg & macs3 bdgcmp -m ppois -t ${out_dir}/${dataset}-pseudoreplicate2_treat_pileup.bdg -c ${out_dir}/${dataset}-pseudoreplicate2_control_lambda.bdg -o ${out_dir}/${dataset}-pseudoreplicate2_ppois.bdg & macs3 bdgcmp -m ppois -t ${out_dir}/${dataset}-pseudoreplicateT_treat_pileup.bdg -c ${out_dir}/${dataset}-pseudoreplicateT_control_lambda.bdg -o ${out_dir}/${dataset}-pseudoreplicateT_ppois.bdg
	echo -e "\t- ${dataset} combining p-value bedgraphs"
	macs3 cmbreps -m fisher -i ${out_dir}/${dataset}-pseudoreplicate1_ppois.bdg ${out_dir}/${dataset}-pseudoreplicate2_ppois.bdg ${out_dir}/${dataset}-pseudoreplicateT_ppois.bdg -o ${out_dir}/${dataset}-combined_ppois.bdg
	sort -k1,1 -k2,2n -S 10% ${out_dir}/${dataset}-combined_ppois.bdg | head -n -1 > ${out_dir}/${dataset}-combined_ppois_sorted.bdg
	echo -e "\t- ${dataset} converting p-value bedgraphs to bigwigs"
	bedGraphToBigWig ${out_dir}/${dataset}-combined_ppois_sorted.bdg ${chr_order} ${out_dir}/${dataset}-pval.bw

	echo "${dataset} completed!"

}
export -f callpeak

datadir="${1}"
chr_order_file="${2}"
blacklist_file="${3}"
parallel="${4}"

frags_dir="${datadir}/pseudobulked_fragments"

ls ${frags_dir} | grep sorted | cut -d "-" -f 1 | sort | uniq | parallel --linebuffer -j ${parallel} callpeak {} ${datadir} ${chr_order_file} ${blacklist_file}
