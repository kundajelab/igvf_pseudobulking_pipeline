#!/bin/bash

# Get current script dir
scriptdir="$(dirname "$(realpath $0)")"

annotations="${1}"
barcode_column="${2}"
basedir="${3}"
outloc="${4}"

# decipher
echo -e "\t- deciphering..."
python ${scriptdir}/decipher_lanemap.py -a "${annotations}" -bc "${barcode_column}" -fd "${basedir}/raw_fragments" -o "${outloc}"
echo -e "\tcomplete!"
