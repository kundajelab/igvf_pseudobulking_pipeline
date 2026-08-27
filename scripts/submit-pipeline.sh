#!/usr/bin/env bash
set -euo pipefail

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
repo_dir=$(dirname "$script_dir")
pushd &> /dev/null "$repo_dir"

# The head process only orchestrates -- it is a mostly idle JVM -- so it wants a long walltime on a
# non-preemptible partition rather than much cpu or memory. This is NOT the queue the tasks run on;
# run-pipeline.sh picks that (owners by default).
# NOTE: "dev" cannot be used here. It rejects batch scripts and caps walltime at 2 hours.
partition="normal"
time_limit="1-00:00:00"
cpus=2
mem="8G"
log_dir="$HOME/logs"
dry_run=false
queue="owners"
profile="$(scripts/get-default-profile.sh)"
workspace="$(scripts/get-default-workspace.sh)"
mode="prod"

function usage {
    cat << EOF
Usage: $(basename "$0") [ARGS] metadata [NEXTFLOW_ARGS]

Submit the pipeline to SLURM as a batch job, so the nextflow head process runs on a compute node
instead of a login node. Returns as soon as the job is queued.

* If metadata is specified by alias or accession ID, download it and infer principal analysis accession if needed.
* If metadata is specified by file path, accession must be specified.
* If no metadata is specified, run small(ish) test set.
* Any additional items after metadata are passed directly to nextflow

ARGS:
    -h|--help: Show this message and exit.
    -P|--partition: Partition to run the nextflow head process on. Default: $partition
    -T|--time: Walltime for the head process, which caps the whole pipeline. Default: $time_limit
    -C|--cpus: CPUs for the head process. Default: $cpus
    -M|--mem: Memory for the head process. Default: $mem
    -l|--log-dir: Where to write the head process log. Default: $log_dir
    -n|--dry-run: Print the sbatch script instead of submitting it.
    -p|--profile: Use comma-separated nextflow profiles. Defaults to inferred from environment ($profile)
    -q|--queue: If running via SLURM, use this queue. Default: $queue
    -w|--workspace: Where to output files. Defaults to inferred from environment ($workspace)
    -a|--principal-analysis: Specify the accession of the principal analysis set
    -m|--mode: Specify IGVF server: "prod", "staging", or "sandbox" ($mode)
EOF
}

# Collect anything we don't consume ourselves and hand it to run-pipeline.sh. The short flags here
# are deliberately upper case so they cannot collide with run-pipeline.sh's -p/-w/-q/-a.
run_args=()
principal_analysis=""
while [[ "$#" -ge 1 ]]; do
    case "$1" in
        "-h" | "--help")
            usage
            exit 0
            ;;
        "-P" | "--partition")
            partition="$2"
            shift 2
            ;;
        "-T" | "--time")
            time_limit="$2"
            shift 2
            ;;
        "-C" | "--cpus")
            cpus="$2"
            shift 2
            ;;
        "-M" | "--mem")
            mem="$2"
            shift 2
            ;;
        "-l" | "--log-dir")
            log_dir="$2"
            shift 2
            ;;
        "--dry-run")
            dry_run="true"
            shift 1
            ;;
        "--no-dry-run")
            dry_run="false"
            shift 1
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
        --)
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

metadata="${1:-"$repo_dir/test_metadata.tsv"}"
shift 1
nextflow_args="${*}"
if [[ "$metadata" =~ \.tsv(\.gz)?$ ]]; then
    metadata_file="$metadata"
    # ensure we have the principal analysis accession
    if [[ -z "$principal_analysis" ]]; then
        if [[ "$metadata_file" =~ test_metadata\.tsv$ ]]; then
            principal_analysis=IGVFDS5417HJRJ,IGVFDS6430MYNQ
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
fi

run_args=("--profile" "$profile" "--workspace" "$workspace" "--queue" "$queue" "--principal-analysis" "$principal_analysis" "--mode" "$mode" "--" "$metadata" "${nextflow_args[*]}")
run_args_quoted=$(printf '%q ' "${run_args[@]}")

job_name=igvf-pseudobulk/$principal_analysis
log_folder="$log_dir/$job_name"
mkdir -p "$log_folder"

sbatch_script=$(cat << EOF
#!/usr/bin/env bash
#SBATCH --job-name=$job_name
#SBATCH --output="$log_folder/%j.out"
#SBATCH --partition=$partition
#SBATCH --time=$time_limit
#SBATCH --cpus-per-task=$cpus
#SBATCH --mem=$mem
#SBATCH --chdir="$repo_dir"
# Requeueing would silently restart the pipeline from scratch, so leave restarts to a manual
# resubmission, which resumes.
#SBATCH --no-requeue
# Ask for SIGINT two minutes before the walltime so nextflow can cancel its outstanding task jobs
# rather than being SIGKILLed and leaving orphans in the queue.
#SBATCH --signal=B:INT@120
set -euo pipefail

echo "head process: job \$SLURM_JOB_ID on \$(hostname), partition \$SLURM_JOB_PARTITION"

# sbatch propagates this job's environment to the task jobs nextflow submits (process.clusterOptions
# sets --export=ALL), which would make each task believe it is running inside this allocation and
# inherit this job's cpu and memory limits. Drop those, but keep SLURM_CONF: the slurm commands on
# the compute nodes need it to reach the controller.
unset \$(env | grep -o '^SLURM_[^=]*' | grep -v '^SLURM_CONF\$' | tr '\n' ' ') || true

# The log is a file rather than a terminal, so nextflow's live-redraw output would land as thousands
# of ANSI escapes instead of readable progress lines.
export NXF_ANSI_LOG=false
export NXF_OPTS='-Xms512m -Xmx6g'

pixi run pipeline $run_args_quoted
EOF
)

if [[ "$dry_run" == "true" ]]; then
    echo "$sbatch_script"
    exit 0
fi

job_id=$(printf '%s\n' "$sbatch_script" | sbatch --parsable)

cat << EOF
submitted nextflow head process as job $job_id
   partition: $partition   walltime: $time_limit   cpus: $cpus   mem: $mem
   log:       $log_folder/$job_id.out
   status:    scripts/status-pipeline.sh $job_id
EOF
