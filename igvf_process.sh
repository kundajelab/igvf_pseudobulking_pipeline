#!/bin/bash

# Get current script dir
scriptdir="$(dirname "$(realpath $0)")"
echo ${scriptdir}

pipeline="${1}"
shift

case "${pipeline}" in
	download)
		echo "Running download pipeline..."
		bash "${scriptdir}/core/0_download_data/PIPELINE.sh" "$@"
		;;
	decipher)
		echo "Decipher annotation lanes..."
		bash "${scriptdir}/core/decipher/PIPELINE.sh" "$@"
		;;
	pseudobulk)
		echo "Running pseudobulk pipeline..."
		bash "${scriptdir}/core/1_pseudobulking/PIPELINE.sh" "$@"
		;;
	*)
		echo "Invalid pipeline: ${pipeline}"
		echo -e "Usage:\nbash igvf_process.sh\n\tdownload \${accession_list} \${access_key} \${secret_key} \${workspace}\n\tdecipher \${annotation_file} \${barcode_lane_column} \${workspace}\n\tpseudobulk \${workspace} \${metadata_path} \${chr_order_file} \${blacklist_file} \${numcpus}"
		exit 1
		;;
esac

