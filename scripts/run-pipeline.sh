#!/usr/bin/env bash
set -euo pipefail

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
repo_dir=$(dirname "$script_dir")
pushd &> /dev/null "$repo_dir"

queue="owners"
profile="$(scripts/get-default-profile.sh)"
workspace="$(scripts/get-default-workspace.sh)"
mode="prod"

function usage {
    cat << EOF
Usage: $0 [ARGS] -- [metadata] [nextflow_args]

Run the pipeline with specified metadata.
* If metadata is specified by alias or accession ID, download it and infer principal analysis accession if needed.
* If metadata is specified by file path, accession must be specified.
* If no metadata is specified, run small(ish) test set.
* Any additional items after metadata are passed directly to nextflow

ARGS:
    -h|--help: Show this message and exit.
    -p|--profile: Use comma-separated nextflow profiles. Defaults to inferred from environment ($profile)
    -q|--queue: If running via SLURM, use this queue. Default: $queue
    -w|--workspace: Where to output files. Defaults to inferred from environment ($workspace)
    -a|--principal-analysis: Specify the accession of the principal analysis set
    -m|--mode: Specify IGVF server: "prod", "staging", or "sandbox" ($mode)
EOF
}

principal_analysis=""
while [[ "$#" -ge 1 ]]; do
    case "$1" in
        "-h" | "--help")
            usage
            exit 0
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
        "-a" | "--principal-analysis")
            principal_analysis="$2"
            shift 2
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

pushd "$repo_dir" &> /dev/null

metadata="${1:-"$repo_dir/test_metadata.tsv"}"
if [[ "$#" -ge 1 ]]; then
    shift 1
fi
nextflow_args="${*}"
if [[ "$metadata" =~ \.tsv(\.gz)?$ ]]; then
    metadata_file="$metadata"
    # ensure we have the principal analysis accession
    if [[ -z "$principal_analysis" ]]; then
        if [[ "$metadata_file" =~ test_metadata\.tsv$ ]]; then
            principal_analysis=IGVFDS5417HJRJ,IGVFDS6430MYNQ
            run_folder="$workspace/${principal_analysis//,/-}"
            mkdir -p "$run_folder"
        else
            1>&2 echo "Must specify metadata accession, or metadata file and principal analysis accession"
            exit 1
        fi
    fi
else
    # ensure we have the principal analysis accession
    if [[ -z "$principal_analysis" ]]; then
        metadata_record=$(pixi run --manifest-path igvf_portal igvf-portal lookup-record "$metadata" --igvf-mode "$mode")
        content_type=$(pixi run yq -r '.content_type' <<< "$metadata_record")
        if [[ "$content_type" != "cell annotations" ]]; then
            1>&2 echo "Got wrong content type ($content_type) for $metadata".
            exit 1
        fi
        principal_analysis=$(pixi run yq -r '.file_set.accession' <<< "$metadata_record")
        if [[ "$principal_analysis" == "null" ]]; then
            1>&2 echo "Unable to find principal analysis for $metadata."
            exit 1
        fi
    fi
    run_folder="$workspace/${principal_analysis//,/-}"
    # download the metadata file if it isn't already present
    metadata_file="$run_folder/${metadata}.tsv.gz"
    if [[ ! -f "$metadata_file" ]]; then
        pixi run --manifest-path igvf_portal igvf-portal download-file "$metadata" --output "$metadata_file" --igvf-mode "$mode"
    fi
fi

# move to the output folder for this workflow so that nextflow logs and folder are stored there
pushd "$run_folder" &> /dev/null

nextflow run "$repo_dir/main.nf" \
    --metadata_file "$metadata_file" \
    --principal_analysis "$principal_analysis" \
    --workspace "$workspace" \
    --slurm_queue "$queue" \
    -profile "$profile" \
    --igvf_mode "$mode" \
    "${nextflow_args[@]}"
