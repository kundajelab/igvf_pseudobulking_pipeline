#!/usr/bin/env bash
set -euo pipefail

queue="normal"
time_limit="1-00:00:00"
cpus=2
mem="8G"
log_dir="$HOME/logs"
work_dir="$(pwd)"
job_name=""
requeue_arg="--requeue"

function usage {
    cat << EOF
Usage: $(basename "$0") [ARGS] -- command args

Submit a command to SLURM as a batch job, so that it runs on a compute node instead of a login node.
Returns as soon as the job is queued.

ARGS:
    -h|--help: Show this message and exit.
    -T|--time: Walltime for the head process, which caps the whole pipeline. Default: $time_limit
    -C|--cpus: CPUs for the head process. Default: $cpus
    -M|--mem: Memory for the head process. Default: $mem
    -l|--log-dir: Where to write the head process log. Default: $log_dir
    -q|--queue: SLURM, use this slurm queue. Default: $queue
    --requeue/--no-requeue: Change whether preempted jobs are requeued. Default: $requeue_arg
    -W|--work-dir: Where to execute the command. Defaults to pwd ($work_dir)
    -n|--name: Job name
EOF
}

while [[ "$#" -ge 1 ]]; do
    case "$1" in
        "-h" | "--help")
            usage
            exit 0
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
        "-W" | "--work-dir")
            work_dir="$2"
            shift 2
            ;;
        "-q" | "--queue")
            queue="$2"
            shift 2
            ;;
        "--requeue" | "--no-requeue")
            requeue_arg="$1"
            shift 1
            ;;
        "-n" | "--name")
            job_name="$2"
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

if [[ "$#" -lt 1 ]]; then
    1>&2 echo "No command was specified."
    usage
    exit 1
fi
if [[ -z "$job_name" ]]; then
    job_name=$(IFS=- ; echo "$*")
fi

command_quoted=$(printf '%q ' "$@")
mkdir -p "$log_dir/$job_name"

sbatch_script=$(cat << EOF
#!/usr/bin/env bash
#SBATCH --job-name=$job_name
#SBATCH --output="$log_dir/$job_name/%j.out"
#SBATCH --partition=$queue
#SBATCH --time=$time_limit
#SBATCH --cpus-per-task=$cpus
#SBATCH --mem=$mem
#SBATCH --chdir="$work_dir"
# Set requeue policy
#SBATCH $requeue_arg
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

$command_quoted
EOF
)

job_id=$(printf '%s\n' "$sbatch_script" | sbatch --parsable)

cat << EOF
submitted nextflow head process as job $job_id
   partition: $queue   walltime: $time_limit   cpus: $cpus   mem: $mem
   log:       $log_dir/$job_name/$job_id.out
   status:    pixi run status --watch $job_id
EOF
