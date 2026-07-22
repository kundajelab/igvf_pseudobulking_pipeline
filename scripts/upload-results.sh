#!/usr/bin/env bash
set -euo pipefail

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
repo_dir=$(dirname "$script_dir")
pushd &> /dev/null "$repo_dir"

queue="normal"
profile="$(scripts/get-default-profile.sh)"
workspace="$(scripts/get-default-workspace.sh)"
mode="prod"
time_limit="1-00:00:00"
cpus=2
mem="8G"
log_dir="$HOME/logs"

function usage {
    cat << EOF
Usage: $0 [ARGS] -- [metadata]

For a pipeline that succeeded with --dry-run for uploads to the IGVF portal, actually perform the upload.
Note: this uses scripts/run-in-container.sh which currently only works with apptainer, so it must be run
on sherlock.

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

metadata="$1"
if [[ "$metadata" =~ \.tsv(\.gz)?$ ]]; then
    metadata_file="$metadata"
    # ensure we have the principal analysis accession
    if [[ -z "$principal_analysis" ]]; then
        1>&2 echo "Must specify metadata accession, or metadata file and principal analysis accession"
        exit 1
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

fi
if [[ ! -f "$metadata_file" ]]; then
    1>&2 echo "Unable to find metadata file '$metadata_file'."
    exit 1
fi

output_folder="$run_folder/output"
dry_run_upload_script="$output_folder/upload.sh"
earnest_upload_script="$output_folder/upload-no-dry-run.sh"
sed 's/^dry_run_arg=".*"$/dry_run_arg=""/' "$dry_run_upload_script" > "$earnest_upload_script"

# Need to run upload in a non-preemptible queue, not the login node
job_name=igvf-upload/$principal_analysis
log_folder="$log_dir/$job_name"
mkdir -p "$log_folder"
sbatch_script=$(cat << EOF
#!/usr/bin/env bash
#SBATCH --job-name=$job_name
#SBATCH --output=$log_folder/%j.out
#SBATCH --partition=$queue
#SBATCH --time=$time_limit
#SBATCH --cpus-per-task=$cpus
#SBATCH --mem=$mem
#SBATCH --chdir=$repo_dir
set -euo pipefail

echo "head process: job \$SLURM_JOB_ID on \$(hostname), partition \$SLURM_JOB_PARTITION"

pixi run-in-container --project igvf_portal "$earnest_upload_script"
EOF
)

job_id=$(printf '%s\n' "$sbatch_script" | sbatch --parsable)

cat << EOF
submitted nextflow head process as job $job_id
   partition: $queue   walltime: $time_limit   cpus: $cpus   mem: $mem
   log:       $log_folder/$job_id.out
   status:    scripts/status-pipeline.sh $job_id
EOF
