#!/usr/bin/env bash
set -euo pipefail

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
repo_dir=$(dirname "$script_dir")

num_simultaneous=2
partition="normal"
time_limit="1-00:00:00"
cpus=2
mem="8G"
log_dir="$HOME/logs"
queue="owners"
profile="$("$script_dir/get-default-profile.sh")"
workspace="$("$script_dir/get-default-workspace.sh")"
mode="prod"
dry_run_arg="--no-dry-run"

function usage {
    cat << EOF
Usage: $(basename "$0") [SUBMIT_ARGS] metadata_id_list [NEXTFLOW_ARGS]

Run the pipeline once per accession as a SLURM array job, so each nextflow head process runs on a
compute node instead of a login node, and at most --num-simultaneous of them run at a time.

metadata_id_list is a list of accession IDs for cell annotations files, preferrably selected from rows of
the output of igvf-portal make-pseudobulk-tracker that are deemed runable (yellow on the google sheet).

SUBMIT_ARGS:
    -h|--help: Show this message and exit.
    -s|--num-simultaneous: Specify number of simultaneous pipelines to run ($num_simultaneous)
    -P|--partition: Partition to run the nextflow head processes on. Default: $partition
    -T|--time: Walltime per pipeline, which caps each whole pipeline. Default: $time_limit
    -C|--cpus: CPUs per head process. Default: $cpus
    -M|--mem: Memory per head process. Default: $mem
    -l|--log-dir: Where to write the head process logs. Default: $log_dir
    -p|--profile: Use comma-separated nextflow profiles. Defaults to inferred from environment ($profile)
    -q|--queue: SLURM, use this slurm queue. Default: $queue
    -w|--workspace: Where to output files. Defaults to inferred from environment ($workspace)
    -m|--mode: Specify IGVF server: "prod", "staging", or "sandbox" ($mode)
    --dry-run/--no-dry-run: Turn on/off igvf_dry_run. With dry-run *on* no changes to the IGVF portal are made.
      With dry-run *off* records are posted/patched and files are uploaded. Default: $dry_run_arg

NEXTFLOW_ARGS are passed through to nextflow by run-pipeline.sh, which each array task runs.

EOF
}

# Collect anything we don't consume ourselves and hand it to run-pipeline.sh. The short flags here
# are deliberately upper case so they cannot collide with run-pipeline.sh's -p/-w/-q/-a.
while [[ "$#" -ge 1 ]]; do
    case "$1" in
        "-h" | "--help")
            usage
            exit 0
            ;;
        "-s" | "--num-simultaneous")
            num_simultaneous="$2"
            shift 2
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
        "--dry-run"|"--no-dry-run")
            dry_run_arg=$1
            shift 1
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

if [[ $# -lt 1 ]]; then
    1>&2 echo "metadata_id_list TSV must be supplied."
    1>&2 usage
    exit 1
fi
metadata_id_list="$1"
if [[ ! -f "$metadata_id_list" ]]; then
    1>&2 echo "metadata_id_list must be a file with a list of metadata accessions to process."
    1>&2 usage
    exit 1
fi
if ! grep -q IGVF_API_KEY <<< "$(env)" || ! grep -q IGVF_SECRET_KEY <<< "$(env)"; then
    1>&2 echo "IGVF_API_KEY and IGVF_SECRET_KEY must be set."
    exit 1
fi
# Validate before queueing anything. A list of the wrong accession type otherwise trickles through
# the whole array at the throttle rate, failing every task with the same error. Flag IGVF accessions
# that are not files, rather than requiring IGVFFI, so aliases (not IGVF-prefixed) still pass.
wrong_type=$(grep -nE '^IGVF' "$metadata_id_list" | grep -vE '^[0-9]+:IGVFFI' || true)
if [[ -n "$wrong_type" ]]; then
    wrong_count=$(printf '%s\n' "$wrong_type" | wc -l | tr -d ' ')
    {
        echo "metadata_id_list must hold cell annotations file accessions (IGVFFI...), one per line."
        echo "$wrong_count line(s) name a different IGVF object type:"
        printf '%s\n' "$wrong_type" | head -5 | sed 's/^\([0-9]*\):/  line \1: /'
        if [[ "$wrong_count" -gt 5 ]]; then
            echo "  ... and $((wrong_count - 5)) more"
        fi
        echo "Use the annotation file accession column of the make-pseudobulk-tracker output, not"
        echo "the principal analysis set accession column."
    } 1>&2
    exit 1
fi
num_jobs=$(wc -l < "$metadata_id_list" | tr -d ' ')
shift 1
nextflow_args=("${@}")

# Only the flags run-pipeline.sh accepts. The SLURM resource flags that submit-pipeline.sh would
# have taken are set as this array's own #SBATCH directives instead, because each array task now is
# the nextflow head process rather than something that submits one.
run_args=(\
    "--profile" "$profile" "--workspace" "$workspace" \
    "--queue" "$queue" "--mode" "$mode" "$dry_run_arg"
)
run_args_quoted_str=$(printf '%q ' "${run_args[@]}")
# printf runs its format once even with no arguments, so an empty nextflow_args would still yield
# "''" here, which would make the array in the job script hold one empty string instead of nothing.
nextflow_args_quoted_str=""
if [[ ${#nextflow_args[@]} -gt 0 ]]; then
    nextflow_args_quoted_str=$(printf '%q ' "${nextflow_args[@]}")
fi

job_name=trickle-submit-pipelines
log_folder="$log_dir/$job_name"
mkdir -p "$log_folder"

# Each array task runs one pipeline to completion, so the task gets the head process resources and
# %num_simultaneous caps how many pipelines run at once. An earlier version had the task submit a
# separate job and exit after ~30s, which made the throttle limit the submission rate instead: the
# pipelines it launched outlived the task and accumulated without bound.
# make header setting sbatch job parameters and variables
sbatch_header=$(cat << EOF
#!/usr/bin/env bash
#SBATCH --job-name=$job_name
#SBATCH --array=1-${num_jobs}%${num_simultaneous}
#SBATCH --output="$log_folder/%A_%a.out"
#SBATCH --partition="$partition"
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

metadata_id_list="$metadata_id_list"
run_args=($run_args_quoted_str)
nextflow_args=($nextflow_args_quoted_str)
EOF
)

# Store the task script without expanding variables, so the SLURM_* references below resolve inside
# the running job rather than here.
sbatch_body=$(cat << 'EOF'
echo "head process: job $SLURM_JOB_ID on $(hostname), partition $SLURM_JOB_PARTITION"

# Pick this task's accession before the unset below clears every SLURM_* variable.
metadata=$(tail -n "+$SLURM_ARRAY_TASK_ID" "$metadata_id_list" | head -n1)
1>&2 echo "Running pipeline for metadata=$metadata"

# sbatch propagates this job's environment to the task jobs nextflow submits (process.clusterOptions
# sets --export=ALL), which would make each task believe it is running inside this allocation and
# inherit this job's cpu and memory limits. Drop those, but keep SLURM_CONF: the slurm commands on
# the compute nodes need it to reach the controller.
unset $(env | grep -o '^SLURM_[^=]*' | grep -v '^SLURM_CONF$' | tr '\n' ' ') || true

# The log is a file rather than a terminal, so nextflow's live-redraw output would land as thousands
# of ANSI escapes instead of readable progress lines.
export NXF_ANSI_LOG=false
export NXF_OPTS='-Xms512m -Xmx6g'

# Sherlock runs bash 4.2, where `set -u` treats an empty array as unset, so a bare
# "${nextflow_args[@]}" aborts when no nextflow args were given. The ${x[@]+...} guard expands to
# nothing in that case and is required by bash < 4.4; the outer expansion must stay unquoted.
pixi run pipeline "${run_args[@]}" -- "$metadata" ${nextflow_args[@]+"${nextflow_args[@]}"}
EOF
)

job_id=$(printf '%s\n' "$sbatch_header" "$sbatch_body" | sbatch --parsable)

cat << EOF
submitted $num_jobs pipelines as array job $job_id (max $num_simultaneous running at a time)
   partition: $partition   walltime: $time_limit   cpus: $cpus   mem: $mem   (per pipeline)
   log:       $log_folder/${job_id}_*.out   (one per pipeline)
   status:    scripts/status-pipeline.sh $job_id
EOF
