#!/bin/bash

# Get current script dir
scriptdir="$(dirname "$(realpath $0)")"

outdir="${1}"
accessions="${2}"
access_key="${3}"
secret_key="${4}"

# 0 - set up workspace
echo -e "\t- setting up workspace..."
bash ${scriptdir}/0_set_up_workspace.sh "${outdir}"

# 1 - download jsons
echo -e "\t- downloading jsons..."
python ${scriptdir}/1_download_jsons.py -a "${accessions}" -k "${access_key}" -s "${secret_key}" -o "${outdir}"

# 2 - download files
echo -e "\t- downloading files..."
python ${scriptdir}/2_download_files.py -k "${access_key}" -s "${secret_key}" -o "${outdir}"

# 3 - rename files
echo -e "\t- renaming files..."
python ${scriptdir}/3_rename_files.py -o "${outdir}"

# 4 cleanup workspace
echo -e "\t- cleaning workspace..."
bash ${scriptdir}/4_cleanup_workspace.sh "${outdir}"

echo -e "\tcomplete!"
