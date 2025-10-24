#!/bin/bash

# Get current script dir
scriptdir="$(dirname "$(realpath $0)")"

basedir="${1}"
metadata_loc="${2}"
chr_order_file="${3}"
blacklist_file="${4}"
parallel="${5}"

# 0 - set up workspace
echo -e "\t- setting up workspace..."
bash ${scriptdir}/0_set_up_workspace.sh "${basedir}"
# 1 - split fragments
echo -e "\t- splitting fragments..."
bash ${scriptdir}/1_split_fragments.sh "${basedir}" "${metadata_loc}" "${chr_order_file}" "${parallel}"
# 2 - catsort
echo -e "\t- concatenating+sorting..."
bash ${scriptdir}/2_catsort.sh "${basedir}" "${parallel}"
# 3 - peak calling
echo -e "\t- calling peaks..."
bash ${scriptdir}/3_call_peaks.sh "${basedir}" "${chr_order_file}" "${blacklist_file}" "${parallel}"
# 4 - rna pseudobulking
echo -e "\t- rna pseudobulking..."
python ${scriptdir}/4_rna_pseudobulking.py -d "${basedir}" -m "${metadata_loc}"
# 5 - rename files
echo -e "\t- renaming files..."
bash ${scriptdir}/5_rename_files.sh ${basedir}
# 6 - clean up workspace
echo -e "\t- cleaning up workspace..."
bash ${scriptdir}/6_cleanup_workspace.sh ${basedir}

echo -e "\tcomplete!"
