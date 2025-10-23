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
	pseudobulk)
		echo "Running pseudobulk pipeline..."
		bash "${scriptdir}/core/1_pseudobulking/PIPELINE.sh" "$@"
		;;
	*)
		echo "Invalid pipeline: ${pipeline}"
		exit 1
		;;
esac

