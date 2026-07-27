#!/usr/bin/env bash
set -euo pipefail

log_dir="$HOME/logs"
task_queue="owners"
lines=25
follow=false
job_id=""

function usage {
    cat << EOF
Usage: $0 [ARGS] [job_id]

Show the status of a pipeline submitted with scripts/submit-pipeline.sh: the state of the nextflow
head process, a breakdown of the task jobs it has in the queue, and the tail of its log.

If no job id is given, use the most recent igvf-pseudobulk job, preferring one that is still queued
or running.

ARGS:
    -h|--help: Show this message and exit.
    -q|--queue: Queue the tasks run on, to summarize. Default: $task_queue
    -l|--log-dir: Where submit-pipeline.sh wrote the head process log. Default: $log_dir
    -N|--lines: How many lines of the head process log to show. Default: $lines
    -f|--follow: Follow the log instead of printing the last lines.
EOF
}

while [[ "$#" -ge 1 ]]; do
    case "$1" in
        "-h" | "--help")
            usage
            exit 0
            ;;
        "-q" | "--queue")
            task_queue="$2"
            shift 2
            ;;
        "-l" | "--log-dir")
            log_dir="$2"
            shift 2
            ;;
        "-N" | "--lines")
            lines="$2"
            shift 2
            ;;
        "-f" | "--follow")
            follow=true
            shift 1
            ;;
        *)
            job_id="$1"
            shift 1
            ;;
    esac
done

if [[ -z "$job_id" ]]; then
    # prefer a live job, otherwise fall back to the most recent one in the accounting database
    job_id=$(
        squeue -h -u "$USER" -o "%i %j" \
            | awk '$2 ~ /^igvf-pseudobulk/ {print $1}' \
            | tail -n1
    )
fi
if [[ -z "$job_id" ]]; then
    job_id=$(
        sacct -u "$USER" -S "$(date -d '14 days ago' +%F)" -X -n -P -o JobID,JobName \
            | awk -F'|' '$2 ~ /^igvf-pseudobulk/ {print $1}' \
            | tail -n1
    )
fi
if [[ -z "$job_id" ]]; then
    1>&2 echo "No igvf-pseudobulk job found. Submit one with scripts/submit-pipeline.sh"
    exit 1
fi

# the label in the log name is only known to submit-pipeline.sh, so match on the job id
log_files=("$log_dir"/igvf-pseudobulk-*-"$job_id".out)
log_file=""
if [[ -e "${log_files[0]}" ]]; then
    # NOTE: negative array indices need bash 4.3, and the login nodes have 4.2
    log_file="${log_files[$(( ${#log_files[@]} - 1 ))]}"
fi

# NOTE: for a job that has left the queue, squeue prints "slurm_load_jobs error: Invalid job id
# specified" on stdout and still exits 0, so neither redirecting stderr nor checking the exit status
# detects it. Keep only output that looks like a job state.
head_state=$(squeue -h -j "$job_id" -o "%T" 2>/dev/null | grep -E '^[A-Z_]+$' || true)

echo "=== head process (job $job_id)"
if [[ -n "$head_state" ]]; then
    squeue -j "$job_id" -o "%.12i %.26j %.10P %.9T %.11M %.11l %.5C %.7m %R"
else
    echo "not in the queue -- final state from accounting:"
    sacct -j "$job_id" -X -o "JobID%14,JobName%26,State%22,ExitCode,Elapsed,Timelimit"
fi

# The workflow.onComplete handler in nextflow.config prints a banner as the last thing in the log, so
# report it here rather than making the outcome something you have to spot in the log tail. A
# shutdown WARN can follow the banner, so match on the banner itself instead of reading the end of
# the file.
echo
echo "=== outcome"
if [[ -z "$log_file" ]]; then
    echo "   (no log yet)"
elif grep -q '^PIPELINE ' "$log_file"; then
    sed -n '/^PIPELINE /,/^====*$/{/^====*$/d; p;}' "$log_file" | sed 's/^/   /'
else
    echo "   still running: nextflow has not printed a completion banner yet"
fi

echo
echo "=== task jobs on $task_queue"
task_summary=$(
    squeue -h -u "$USER" -p "$task_queue" -o "%T|%j" \
        | awk -F'|' '{
            name = $2
            sub(/_\([0-9]+\)$/, "", name)
            sub(/^nf-/, "", name)
            print $1, name
        }' \
        | sort \
        | uniq -c \
        | sort -rn
)
if [[ -n "$task_summary" ]]; then
    printf '%s\n' "$task_summary"
    printf '%7d TOTAL\n' "$(squeue -h -u "$USER" -p "$task_queue" -o "%i" | wc -l)"
else
    echo "   (none queued or running)"
fi

echo
if [[ -z "$log_file" ]]; then
    echo "=== head process log: not found under $log_dir for job $job_id"
    exit 0
fi

if [[ "$follow" == true ]]; then
    echo "=== following $log_file"
    tail -f "$log_file"
else
    echo "=== last $lines lines of $log_file"
    tail -n "$lines" "$log_file"
fi
