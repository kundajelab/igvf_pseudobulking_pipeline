#!/usr/bin/env bash
set -euo pipefail

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
repo_dir=$(dirname "$script_dir")

queue="normal"
profile="$(scripts/get-default-profile.sh)"
workspace="$(scripts/get-default-workspace.sh)"
principal_analysis=""
args=()
while [[ "$#" -ge 1 ]]; do
    case "$1" in
        "-p" | "--profile")
            profile="$2"
            shift 2
            ;;
        "-w" | "--workspace")
            workspace="$2"
            shift 2
            ;;
        "-q" | "--queue")
            queue="$2"
            shift 2
            ;;
        "-a" | "--principal-analysis")
            principal_analysis="$2"
            shift 2
            ;;
        "--")
            break
            shift 1
            ;;
        *)
            break
            ;;
    esac
done

metadata_file="$1"
if [[ -z "$principal_analysis" ]]; then    
    principal_analysis=$(uv run --project igvf_portal igvf-portal infer-principal-analysis "$metadata_file")
fi


nextflow run . \
    --metadata_file "$metadata_file" \
    --principal_analysis "$principal_analysis" \
    --workspace "$workspace" \
    --slurm_queue "$queue" \
    -profile "$profile"
