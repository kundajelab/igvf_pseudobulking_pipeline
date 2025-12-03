#!/bin/bash

# Get current script dir
scriptdir="$(dirname "$(realpath $0)")"

basedir="${1}"
metadata_loc="${2}"
chr_sizes="${3}"
blacklist_file="${4}"
tss_file="${5}"
gene_info="${6}"
ncpus="${7}"

echo "basedir: ${basedir}"
echo "metadata: ${metadata_loc}"
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
bash ${scriptdir}/1_split_fragments.sh "${basedir}" "${metadata_loc}" "${chr_sizes}" "${tss_file}" "${ncpus}"
# 2 - catsort
echo -e "\t- concatenating+sorting..."
bash ${scriptdir}/2_catsort.sh "${basedir}" "${ncpus}"
# 3 - peak calling
echo -e "\t- calling peaks..."
# bash ${scriptdir}/3_call_peaks.sh "${basedir}" "${chr_order_file}" "${blacklist_file}" "${ncpus}"
# 4 - rna pseudobulking
echo -e "\t- rna pseudobulking..."
# python ${scriptdir}/4_rna_pseudobulking.py -d "${basedir}" -m "${metadata_loc}"
# 5 - rename files
echo -e "\t- renaming files..."
# bash ${scriptdir}/5_rename_files.sh ${basedir}
# 6 - clean up workspace
echo -e "\t- cleaning up workspace..."
# bash ${scriptdir}/6_cleanup_workspace.sh ${basedir}

echo -e "\tcomplete!"
