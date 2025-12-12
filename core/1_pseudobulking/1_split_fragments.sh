#!/bin/bash

datadir="${1}"
metadata_loc="${2}"
at_annotation_level="${3}"
chr_order_file="${4}"
tss_file="${5}"
parallel="${6}"

raw_frags_dir="${datadir}/raw_fragments"

# Get current script dir
scriptdir="$(dirname "$(realpath $0)")"

ls ${raw_frags_dir}/*.bed.gz | grep -oP '(?<=/)\w+(?=\.)' | xargs -I {} -P ${parallel} python ${scriptdir}/1_split_fragments.py -d ${datadir} -f {} -m "${metadata_loc}" -a "${at_annotation_level}" -c "${chr_order_file}" -t "${tss_file}"
