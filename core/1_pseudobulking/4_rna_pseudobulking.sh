#!/bin/bash

datadir="${1}"
metadata_loc="${2}"

raw_rna_dir="${datadir}/raw_rna"

# Get current script dir
scriptdir="$(dirname "$(realpath $0)")"

ls ${raw_rna_dir}/*.h5ad | grep -oP '(?<=/)\w+(?=\.)' | xargs -I {} python ${scriptdir}/4_rna_pseudobulking.py -d ${datadir} -a {} -m "${metadata_loc}"

