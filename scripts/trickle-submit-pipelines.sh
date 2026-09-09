#!/usr/bin/env bash
set -euo pipefail

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
repo_dir=$(dirname "$script_dir")

num_simultaneous=2
sleep_time=60
partition="normal"
time_limit="1-00:00:00"
cpus=2
mem="8G"
log_dir="$HOME/logs"
dry_run_arg="--no-dry-run"
queue="owners"
profile="$("$script_dir/get-default-profile.sh")"
workspace="$("$script_dir/get-default-workspace.sh")"
mode="prod"

function usage {
    cat << EOF
Usage: $(basename "$0") [SUBMIT_ARGS] pseudobulk_status [NEXTFLOW_ARGS]

Repeatedly submit the pipeline to SLURM as a batch job, so the nextflow head process runs on a compute node
instead of a login node.

pseudobulk_status is a TSV output by igvf-portal make-pseudobulk-tracker, describing each job to submit

SUBMIT_ARGS:
    -h|--help: Show this message and exit.
    -s|--num-simultaneous: Specify number of simultaneous pipelines to run ($num_simultaneous)
    -S|--sleep-time: Wait this long in between checking pipeline status ($sleep_time)

Additional SUBMIT_ARGS are passed to submit-pipeline.sh.
$("$script_dir/submit-pipeline.sh" --help)

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
        "-S" | "--sleep-time")
            sleep_time="$2"
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
        "--dry-run")
            dry_run_arg="--dry-run"
            shift 1
            ;;
        "--no-dry-run")
            dry_run_arg="--no-dry-run"
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

if [[ $# -lt 1 ]]; then
    1>&2 echo "pseudobulk_status TSV must be supplied."
    1>&2 usage
    exit 1
fi
pseudobulk_status="$1"
shift 1
nextflow_args=("${@}")

submit_args=(\
    "--partition" "$partition" "--time" "$time_limit" "--cpus" "$cpus" "--mem" "$mem" "--log-dir"\
    "$log_dir" "$dry_run_arg" "--profile" "$profile" "--workspace" "$workspace" \
    "--queue" "$queue" "--mode" "$mode"
)
submit_args_quoted_str=$(printf '%q ' "${submit_args[@]}")
# printf runs its format once even with no arguments, so an empty nextflow_args would still yield
# "''" here, which would make the array in the job script hold one empty string instead of nothing.
nextflow_args_quoted_str=""
if [[ ${#nextflow_args[@]} -gt 0 ]]; then
    nextflow_args_quoted_str=$(printf '%q ' "${nextflow_args[@]}")
fi

job_name=trickle-submit-pipelines
log_folder="$log_dir/$job_name"
mkdir -p "$log_folder"

# make header setting sbatch job parameters and variables
sbatch_header=$(cat << EOF
#!/usr/bin/env bash
#SBATCH --job-name=$job_name
#SBATCH --output="$log_folder/%j.out"
#SBATCH --partition=akundaje
#SBATCH --time=7:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH --chdir="$repo_dir"
set -euo pipefail

pseudobulk_status="$pseudobulk_status"
num_simultaneous="$num_simultaneous"
sleep_time="$sleep_time"
submit_args=($submit_args_quoted_str)
nextflow_args=($nextflow_args_quoted_str)
EOF
)

# Store the task script without expanding variables, so that awk script works
sbatch_body=$(cat << 'EOF'
echo "head process: job $SLURM_JOB_ID on $(hostname), partition $SLURM_JOB_PARTITION"

# sbatch propagates this job's environment to the task jobs nextflow submits (process.clusterOptions
# sets --export=ALL), which would make each task believe it is running inside this allocation and
# inherit this job's cpu and memory limits. Drop those, but keep SLURM_CONF: the slurm commands on
# the compute nodes need it to reach the controller.
unset $(env | grep -o '^SLURM_[^=]*' | grep -v '^SLURM_CONF$' | tr '\n' ' ') || true

# The log is a file rather than a terminal, so nextflow's live-redraw output would land as thousands
# of ANSI escapes instead of readable progress lines.
export NXF_ANSI_LOG=false
export NXF_OPTS='-Xms512m -Xmx6g'

function filter_jobs() {
    awk -F '\t' -v OFS='\t' -v RS='\r?\n' \
        '
        NR==1 { for(i=1; i<=NF; ++i){ fields[$i] = i }; next }
        {
            pseudobulking_status=$fields["pseudobulking status"]
            missing_annotations_columns=$fields["missing annotations columns"]
            uniform_pipeline_status=$fields["uniform pipeline status"]
            if (pseudobulking_status=="unattempted" && missing_annotations_columns=="" && uniform_pipeline_status == "completed") {
                print $2, $1
            }
        }
        ' \
        "$pseudobulk_status"
}

function get_num_running() {
    squeue -u "$USER" -O "name:50" | grep -c igvf-pseudobulk
}

#shellcheck disable=SC2218
filter_jobs \
| while read -r metadata principal_analysis; do
    while [[ $(get_num_running) -ge $num_simultaneous ]]; do
        sleep "$sleep_time"
    done
    1>&2 echo "Submitting metadata=$metadata, principal_analysis=$principal_analysis"
    pixi run submit "${submit_args[@]}" --principal-analysis "$principal_analysis" -- "$metadata" "${nextflow_args[@]}"
done
EOF
)

job_id=$(printf '%s\n' "$sbatch_header" "$sbatch_body" | sbatch --parsable)

cat << EOF
submitted nextflow head process as job $job_id
   partition: $partition   walltime: $time_limit   cpus: $cpus   mem: $mem
   log:       $log_folder/$job_id.out
   status:    scripts/status-pipeline.sh $job_id
EOF
