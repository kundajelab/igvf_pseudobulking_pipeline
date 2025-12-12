#!/bin/bash

# Get current script dir
scriptdir="$(dirname "$(realpath $0)")"

basedir="${1}"
metadata_loc="${2}"
at_annotation_level="${3}"
chr_sizes="${4}"
blacklist_file="${5}"
tss_file="${6}"
gene_info="${7}"
ncpus="${8}"

echo "basedir: ${basedir}"
echo "metadata: ${metadata_loc}"
echo "at_annotation_level: ${at_annotation_level}"
echo "chr sizes: ${chr_sizes}"
echo "blacklist_file: ${blacklist_file}"
echo "tss file: ${tss_file}"
echo "gene info: ${gene_info}"
echo "ncpus: ${ncpus}"

# 0 - set up workspace
echo -e "\t- setting up workspace..."
bash ${scriptdir}/0_set_up_workspace.sh "${basedir}"

# 1 - split fragments
echo -e "\t- splitting fragments..."
bash ${scriptdir}/1_split_fragments.sh "${basedir}" "${metadata_loc}" "${at_annotation_level}" "${chr_sizes}" "${tss_file}" "${ncpus}"

# 2 - catsort
echo -e "\t- concatenating+sorting..."
bash ${scriptdir}/2_catsort.sh "${basedir}" "${ncpus}" "${chr_sizes}"

# 3 - peak calling
echo -e "\t- calling peaks..."
bash ${scriptdir}/3_call_peaks.sh "${basedir}" "${chr_sizes}" "${blacklist_file}" "${ncpus}"

# 4 - rna pseudobulking
echo -e "\t- rna pseudobulking..."
python ${scriptdir}/4_rna_pseudobulking.py -d "${basedir}" -m "${metadata_loc}" -a "${at_annotation_level}" -g "${gene_info}"

# 5 - aggregate_pseudobulk_outputs
echo -e "\t- aggregating pseudobulk outputs..."
python ${scriptdir}/5_aggregate_outputs.py -d ${basedir} -m "${metadata_loc}" -a "${at_annotation_level}"

if [ -f "${basedir}/step5_complete.txt" ]; then
    echo -e "\t- cleaning up workspace..."
    bash "${scriptdir}/6_cleanup_workspace.sh" "${basedir}"
    echo -e "\tcomplete!"
else
    echo "ERROR FOUND. NOT CLEANING UP DIRECTORY."
    exit 1
fi