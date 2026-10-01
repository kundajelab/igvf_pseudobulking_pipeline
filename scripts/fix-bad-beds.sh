#!/usr/bin/env bash
set -euo pipefail

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
repo_dir=$(dirname "$script_dir")
pushd &> /dev/null "$repo_dir"

queue="owners"
profile="$(scripts/get-default-profile.sh)"
workspace="$(scripts/get-default-workspace.sh)"
mode="prod"
extra_args=()

function usage {
    cat << EOF
Usage: $0 [ARGS] -- [nextflow_args]

Run the fix_bad_beds.nf pipeline.
* Any additional items are passed directly to nextflow

ARGS:
    -h|--help: Show this message and exit.
    -l|--lab: Fix bad beds for this lab only
    -p|--profile: Use comma-separated nextflow profiles. Defaults to inferred from environment ($profile)
    -q|--queue: If running via SLURM, use this queue. Default: $queue
    -w|--workspace: Where to output files. Defaults to inferred from environment ($workspace)
    -m|--mode: Specify IGVF server: "prod", "staging", or "sandbox" ($mode)
    --max-pseudobulks: Specify the maximum number of pseudobulks to fix. (unlimited)
EOF
}

while [[ "$#" -ge 1 ]]; do
    case "$1" in
        "-h" | "--help")
            usage
            exit 0
            ;;
        "-l" | "--lab")
            extra_args+=("--lab" "$2")
            shift 2
            ;;
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
        "--dry-run")
            extra_args+=("--dry_run" "true")
            shift 1
            ;;
        "--no-dry-run")
            extra_args+=("--dry_run" "false")
            shift 1
            ;;
        "-m" | "--mode")
            mode="$2"
            case "$mode" in
                prod|staging|sandbox)
                    ;;
                *)
                    1>&2 echo "Invalid IGVF portal mode: $mode"
                    exit 1
                    ;;
            esac
            shift 2
            ;;
        --max-pseudobulks)
            extra_args+=(--max_pseudobulks "$2")
            shift 2
            ;;
        "--")
            shift 1
            break
            ;;
        "--"?*)
            1>&2 echo "Unknown argument: $1"
            exit 1
            ;;
        *)
            break
            ;;
    esac
done

run_folder="$workspace/fix_bad_beds"
mkdir -p "$run_folder"
pushd "$run_folder" &> /dev/null

set +x
nextflow \
    -c "$repo_dir/fix_bad_beds.config" \
    run \
    "$repo_dir/fix_bad_beds.nf" \
    -profile "$profile" \
    --workspace "${workspace}" \
    --igvf-mode "$mode" \
    --slurm_queue "${queue}" \
    ${extra_args[@]+"${extra_args[@]}"}
