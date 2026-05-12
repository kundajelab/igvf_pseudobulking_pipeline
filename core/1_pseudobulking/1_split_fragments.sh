#!/bin/bash
set -euo pipefail

datadir="${1}"
metadata_loc="${2}"
chr_order_file="${3}"
tss_file="${4}"
parallel="${5}"

raw_frags_dir="${datadir}/raw_fragments"

# Get current script dir
scriptdir="$(dirname "$(realpath $0)")"

find ${raw_frags_dir} -name "*.bed.gz" | xargs -I {} basename {} .bed.gz \
| xargs -I {} -P ${parallel} python ${scriptdir}/1_split_fragments.py -d ${datadir} -f {} -m "${metadata_loc}" -a "${at_annotation_level}" -c "${chr_order_file}" -t "${tss_file}"

# Validate: every raw fragment file should have a corresponding QC report
failed=0
for frag_file in ${raw_frags_dir}/*.bed.gz; do
    accession=$(basename "${frag_file}" .bed.gz)
    if [ ! -f "${datadir}/atac_qc_reports/${accession}.tsv" ]; then
        echo "ERROR: Step 1 - Missing QC report for ${accession}" >&2
        failed=1
    fi
    if [ ! -f "${datadir}/atac_qc_reports/${accession}_tss_matrix.npy" ]; then
        echo "ERROR: Step 1 - Missing TSS matrix for ${accession}" >&2
        failed=1
    fi
done

if [ ${failed} -eq 0 ]; then
    touch "${datadir}/step1_complete.txt"
    echo "Step 1 (split fragments) completed successfully."
else
    echo "ERROR: Step 1 (split fragments) failed validation." >&2
    exit 1
fi
