#!/usr/bin/env bash

access_key=""
secret_key=""
igvffs_in=""
outdir=""

usage() {
  echo "Usage $0 --access_key <value> --secret_key <value> --igvffs_in <value> --outdir <value>"
  echo "  --access_key: from IGVF portal (see https://data.igvf.org/user-profile/)"
  echo "  --secret_key: from IGVF portal (see https://data.igvf.org/user-profile/)"
  echo "  --igvffs_in: json with accessions for downloading files (see get_igvf_download_accessions.py)"
  echo "  --outdir: directory to store outputs"
}

while [[ $# -gt 0 ]]; do
  case $1 in
    --access_key)
      if [ -z "$2" ] || [[ $2 == --* ]]; then
                echo "Error: --access_key requires a value"
                usage
            fi
            access_key="$2"
            shift 2
            ;;
    --secret_key)
      if [ -z "$2" ] || [[ $2 == --* ]]; then
                echo "Error: --secret_key requires a value"
                usage
            fi
            secret_key="$2"
            shift 2
            ;;
    --igvffs_in)
      if [ -z "$2" ] || [[ $2 == --* ]]; then
                echo "Error: --igvffs requires a value"
                usage
            fi
            igvffs_in="$2"
            shift 2
            ;;
    --outdir)
      if [ -z "$2" ] || [[ $2 == --* ]]; then
                echo "Error: --outdir requires a value"
                usage
            fi
            outdir="$2"
            shift 2
            ;;
    *)
	    echo "Error: invalid usage"
	    usage
	    exit 1
	    ;;
  esac
done

if [ -z "$access_key" ] || [ -z "$secret_key" ] || [ -z "$igvffs_in" ] || [ -z "$outdir" ]; then
    echo "Error: All arguments are required"
    usage
    exit
fi

log="$outdir"/download.log

echo "Collecting IGVFFs from $igvffs_in ..." | tee -a "$log"

# read in values from json
scRNA_igvff=$(cat "$igvffs_in" | jq '."counts_matrix"' | sed -e 's/\"//g')
scATAC_igvff=$(cat "$igvffs_in" | jq '."fragments"' | sed -e 's/\"//g')

echo "scRNA: $scRNA_igvff" | tee -a "$log"
echo "scATAC: $scATAC_igvff" | tee -a "$log"

scATAC_url=https://api.data.igvf.org/tabular-files/"$scATAC_igvff"/@@download/"$scATAC_igvff".bed.gz
echo "Downloading scATAC from $scATAC_url ..." | tee -a "$log"

# download scATAC output
curl -sRL -u "$access_key":"$secret_key" "$scATAC_url" -o "$outdir"/"$scATAC_igvff".bed.gz

scRNA_url=https://api.data.igvf.org/matrix-files/"$scRNA_igvff"/@@download/"$scRNA_igvff".h5ad
echo "Downloading scRNA from $scRNA_url ..." | tee -a "$log"

# download scRNA output
curl -sRL -u "$access_key":"$secret_key" "$scRNA_url" -o "$outdir"/"$scRNA_igvff".h5ad
