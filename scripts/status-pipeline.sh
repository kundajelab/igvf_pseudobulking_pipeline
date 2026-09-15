#!/usr/bin/env bash
set -euo pipefail

log_dir="$HOME/logs"
task_queue="owners"
lines=25
follow=false
watch=false
# The Sherlock login profile caches squeue output for 20s and sacct for 60s, so refreshing more often
# only redraws cached data.
interval=60
verbose=false
max_rows=20
job_id=""

function usage {
    cat << EOF
Usage: $0 [ARGS] [job_id]

Show the status of a pipeline run.

For a job from scripts/submit-pipeline.sh, or a single task given as ARRAYID_TASK: the state of the
nextflow head process, its outcome, a breakdown of the task jobs in the queue, and the tail of its log.

For an array job from scripts/trickle-submit-pipelines.sh: counts per state, then every task that
needs attention (failed, cancelled, timed out) with its accession, the nextflow process that failed
and where in the log, then the running tasks. Succeeded tasks are counted, and listed with -v.

If no job id is given, use the most recent igvf-pseudobulk job or trickle-submit-pipelines array,
preferring one that is still queued or running.

ARGS:
    -h|--help: Show this message and exit.
    -q|--queue: Queue the tasks run on, to summarize. Default: $task_queue
    -l|--log-dir: Where the head process logs were written. Default: $log_dir
    -N|--lines: How many lines of the head process log to show. Default: $lines
    -f|--follow: Follow the log(s) after printing the summary.
    -W|--watch: Refresh the summary until the job or array finishes. For an array, newly failed or
      cancelled tasks are announced, with a terminal bell.
    -i|--interval SECONDS: How often --watch refreshes. Default: $interval. On Sherlock, sacct results
      are cached for 60s, so a shorter interval mostly redraws the same data.
    -v|--verbose: For an array, list every task needing attention rather than the first $max_rows,
      and list the pending and succeeded tasks too.
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
        "-W" | "--watch")
            watch=true
            shift 1
            ;;
        "-i" | "--interval")
            interval="$2"
            shift 2
            ;;
        "-v" | "--verbose")
            verbose=true
            shift 1
            ;;
        "--"?*)
            1>&2 echo "Unknown argument: $1"
            exit 1
            ;;
        *)
            job_id="$1"
            shift 1
            ;;
    esac
done

if [[ ! "$interval" =~ ^[1-9][0-9]*$ ]]; then
    1>&2 echo "--interval needs a whole number of seconds, got: $interval"
    exit 1
fi
if [[ "$watch" == true && "$follow" == true ]]; then
    1>&2 echo "--watch and --follow cannot be combined"
    exit 1
fi

job_pattern='^(igvf-pseudobulk|trickle-submit-pipelines)'
if [[ -z "$job_id" ]]; then
    # Prefer a live job. %F is the array job id (the plain job id for a non-array job), so every task
    # of an array collapses to the id this script wants. Job ids increase, so the largest is newest.
    job_id=$(
        squeue -h -u "$USER" -o "%F %j" \
            | awk -v p="$job_pattern" '$2 ~ p {print $1}' \
            | sort -n \
            | tail -n1
    )
fi
if [[ -z "$job_id" ]]; then
    job_id=$(
        sacct -u "$USER" -S "$(date -d '14 days ago' +%F)" -X -n -P -o JobID,JobName \
            | awk -F'|' -v p="$job_pattern" '$2 ~ p {id = $1; sub(/_.*/, "", id); print id}' \
            | sort -n \
            | tail -n1
    )
fi
if [[ -z "$job_id" ]]; then
    1>&2 echo "No igvf-pseudobulk job or trickle-submit-pipelines array found in the last 14 days."
    exit 1
fi

# the label in the log name is only known to the submitting script, so match on the job id
job_name=$(sacct -j "$job_id" --format=JobName%100 -n | awk 'NR==1{$1=$1;print}')
if [[ -z "$job_name" ]]; then
    1>&2 echo "No accounting record for job $job_id."
    exit 1
fi
log_folder="$log_dir/$job_name"

# A task id (ARRAYID_TASK) has its own log, so only a bare array id gets the array view. awk reads all
# of sacct's output, where grep -q would exit early and fail the pipeline through SIGPIPE.
is_array=false
if [[ "$job_id" != *_* ]] \
    && sacct -j "$job_id" -X -n -P -o JobID \
        | awk -v id="${job_id}_" 'index($0, id) == 1 {found = 1} END {exit !found}'; then
    is_array=true
fi

function clear_screen {
    if [[ -t 1 ]]; then
        printf '\033[H\033[2J'
    else
        echo "----- $(date '+%F %T') -----"
    fi
}

function job_state {
    # NOTE: for a job that has left the queue, squeue reports "slurm_load_jobs error: Invalid job id
    # specified". On Sherlock the login profile wraps squeue in a cache that merges stderr into
    # stdout and loses the exit status, so neither redirecting stderr nor checking the exit status
    # detects it. Keep only output that looks like a job state.
    squeue -h -j "$1" -o "%T" 2> /dev/null | grep -E '^[A-Z_]+$' || true
}

function render_queue {
    # Task jobs are named nf-PROCESS_(N); count them by state and process across all pipelines.
    local task_summary
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
    echo "=== task jobs on $task_queue (all of your pipelines)"
    if [[ -n "$task_summary" ]]; then
        printf '%s\n' "$task_summary"
        printf '%7d TOTAL\n' "$(squeue -h -u "$USER" -p "$task_queue" -o "%i" | wc -l)"
    else
        echo "   (none queued or running)"
    fi
}

function render_single {
    local log_file="$log_folder/$job_id.out"
    [[ -f "$log_file" ]] || log_file=""

    echo "=== head process (job $job_id)"
    if [[ -n "$(job_state "$job_id")" ]]; then
        squeue -j "$job_id" -o "%.12i %.26j %.10P %.9T %.11M %.11l %.5C %.7m %R"
    else
        echo "not in the queue -- final state from accounting:"
        sacct -j "$job_id" -X -o "JobID%14,JobName%26,State%22,ExitCode,Elapsed,Timelimit"
    fi

    # The workflow.onComplete handler in main.nf prints a banner as the last thing in the log, so
    # report it here rather than making the outcome something you have to spot in the log tail. A
    # shutdown WARN can follow the banner, so match on the banner itself instead of reading the end
    # of the file.
    echo
    echo "=== outcome"
    if [[ -z "$log_file" ]]; then
        echo "   (no log yet)"
    elif grep -q '^PIPELINE ' "$log_file"; then
        sed -n '/^PIPELINE /,/^====*$/{/^====*$/d; p;}' "$log_file" | sed 's/^/   /'
    else
        echo "   still running: nextflow has not printed a completion banner yet"
    fi
    if [[ -n "$log_file" ]]; then
        local first_error
        first_error=$(grep -n -m1 -E 'Error executing process > |^ERROR ~ ' "$log_file" || true)
        if [[ -n "$first_error" ]]; then
            echo "   first error, line ${first_error%%:*}: ${first_error#*:}"
        fi
    fi

    echo
    render_queue

    echo
    if [[ -z "$log_file" ]]; then
        echo "=== head process log: not found under $log_dir for job $job_id"
    else
        echo "=== last $lines lines of $log_file"
        tail -n "$lines" "$log_file"
    fi
}

# One pass over every task log, printing one tab-separated line per log:
#   task, log name, accession, banner, task counts, failed process, its exit status, error line,
#   other error message, retries, last process submitted, SLURM cancel reason, last line.
# q is a single quote, which cannot appear inside this single-quoted program.
# shellcheck disable=SC2016
log_awk='
function clean(s) {
    gsub(/[\t\r]/, " ", s)
    sub(/^ +/, "", s)
    sub(/ +$/, "", s)
    return s
}
function flush() {
    if (file == "") return
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", task, base, acc, banner, counts, \
        errproc, errexit, errline, errmsg, retries, lastproc, cancel, lastline
}
FNR == 1 {
    flush()
    file = FILENAME
    base = file
    sub(/.*\//, "", base)
    task = base
    sub(/\.out$/, "", task)
    sub(/^[^_]*_/, "", task)
    acc = banner = counts = errproc = errexit = errline = errmsg = lastproc = cancel = lastline = ""
    retries = 0
    want_exit = 0
}
acc == "" && /^(Running pipeline for|Submitting) metadata=/ {
    acc = $0
    sub(/^[^=]*=/, "", acc)
    acc = clean(acc)
}
/^PIPELINE (SUCCEEDED|FAILED)/ { banner = $2 }
/^ *tasks: +[0-9]+ succeeded/ {
    s = $0
    sub(/^ *tasks: +/, "", s)
    split(s, c, /, */)
    ok = c[1] + 0
    cached = c[2] + 0
    bad = c[3] + 0
    counts = ok " ok" (cached > 0 ? " (" cached " cached)" : "") ", " bad " failed"
}
errproc == "" && index($0, "Error executing process > ") {
    p = substr($0, index($0, "Error executing process > ") + 26)
    gsub(q, "", p)
    errproc = clean(p)
    errline = FNR
    want_exit = 1
}
want_exit && /terminated with an error exit status \(/ {
    x = $0
    sub(/.*exit status \(/, "", x)
    sub(/\).*/, "", x)
    errexit = x
    want_exit = 0
}
errproc == "" && errmsg == "" && /^ERROR ~ / {
    errmsg = clean(substr($0, 9))
    errline = FNR
}
/Execution is retried/ { retries++ }
/(Submitted|Re-submitted) process > / {
    p = $0
    sub(/.*process > /, "", p)
    lastproc = clean(p)
}
/\*\*\* JOB [0-9]+ ON .* CANCELLED AT / {
    x = $0
    sub(/.* CANCELLED AT [^ ]+ /, "", x)
    sub(/ *\*\*\*.*/, "", x)
    cancel = clean(x)
}
NF && !/^[=-]+$/ { lastline = clean($0) }
END { flush() }
'

# Join sacct records (S lines) and squeue tasks (Q lines), both "|"-separated, with the task log
# summaries (L lines, tab-separated) into one tab-separated row per task, or per range of tasks that
# have not started:
#   kind (T task / R range), task or range, category, state, elapsed (or task count for R),
#   accession, detail, log pointer, reason (for grouping identical failures).
# sacct is authoritative for tasks that have started: one that was cancelled, timed out or ran out
# of memory never printed a banner, so the log alone cannot say what happened to it. squeue is
# authoritative for pending tasks while the array is live: the scheduler splits the next tasks due
# to start off the array into their own records, which sacct does not show until they start.
# shellcheck disable=SC2016
join_awk='
# Count the indices in a range spec such as "6-20" or "3,5-7", marking each as seen.
function count_indices(spec,   parts, n, i, j, ab, total) {
    n = split(spec, parts, ",")
    total = 0
    for (i = 1; i <= n; i++) {
        if (split(parts[i], ab, "-") != 2) ab[2] = ab[1]
        for (j = ab[1] + 0; j <= ab[2] + 0; j++) see(j)
        total += ab[2] - ab[1] + 1
    }
    return total
}
function see(index_) {
    seen[index_] = 1
    if (index_ > max_index) max_index = index_
}
function user_name(uid,   cmd, line, f) {
    if (uid in users) return users[uid]
    cmd = "getent passwd " uid
    line = ""
    cmd | getline line
    close(cmd)
    split(line, f, ":")
    users[uid] = (f[1] != "") ? f[1] : "uid " uid
    return users[uid]
}
function exit_note(code) {
    if (code == 137) return " (killed: likely out of memory)"
    if (code == 143) return " (SIGTERM: preempted or cancelled)"
    return ""
}
function mask(s) {
    gsub(/IGVF[A-Z][A-Z][A-Z0-9]+/, "IGVF<ACC>", s)
    return s
}
function pending_detail(why) {
    if (why == "JobArrayTaskLimit") return "waiting for a slot (--num-simultaneous)"
    if (why == "None") return "waiting to be scheduled"
    return "queued: " why
}
function range_spec(first, last) {
    return (first == last) ? first : first "-" last
}
function out(kind, task, cat, state, elapsed, accession, detail, pointer, reason) {
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", kind, task, cat, state, elapsed, accession, \
        detail, pointer, reason
}
{
    tag = substr($0, 1, 1)
    rest = substr($0, 3)
    if (tag == "L") {
        split(rest, f, "\t")
        t = f[1]
        have_log[t] = 1
        base[t] = f[2]; acc[t] = f[3]; banner[t] = f[4]; counts[t] = f[5]
        errproc[t] = f[6]; errexit[t] = f[7]; errline[t] = f[8]; errmsg[t] = f[9]
        retries[t] = f[10]; lastproc[t] = f[11]; cancel[t] = f[12]; lastline[t] = f[13]
    } else if (tag == "S") {
        split(rest, f, "|")
        nrec++
        rid[nrec] = f[1]; rstate[nrec] = f[2]; rexit[nrec] = f[3]; relapsed[nrec] = f[4]
    } else if (tag == "Q") {
        split(rest, f, "|")
        nq++
        qidx[nq] = f[1] + 0; qstate[nq] = f[2]; qelapsed[nq] = f[3]; qreason[nq] = f[4]
    }
}
END {
    live = (nq > 0)
    for (i = 1; i <= nrec; i++) {
        if (index(rid[i], "[")) continue
        t = rid[i]
        sub(/^[^_]*_/, "", t)
        in_sacct[t + 0] = 1
    }
    for (i = 1; i <= nq; i++) {
        k = qidx[i]
        if (qstate[i] == "PENDING") {
            pend[k] = qreason[i]
            see(k)
            if (k > max_pend) max_pend = k
        } else if (!(k in in_sacct)) {
            # started so recently that sacct has no record yet
            nrec++
            rid[nrec] = job "_" k; rstate[nrec] = qstate[i]; rexit[nrec] = "0:0"
            relapsed[nrec] = qelapsed[i]
        }
    }

    for (i = 1; i <= nrec; i++) {
        id = rid[i]
        state = rstate[i]
        who = ""
        if (match(state, / by [0-9]+/)) who = user_name(substr(state, RSTART + 4, RLENGTH - 4))
        sub(/ .*/, "", state)
        # while the array is live, pending tasks come from squeue, which sees all of them
        if (live && state == "PENDING") continue

        if (index(id, "[")) {
            spec = id
            sub(/^[^[]*\[/, "", spec)
            sub(/\].*$/, "", spec)
            sub(/%.*/, "", spec)
            n = count_indices(spec)
            by = (who != "") ? " by " who : ""
            if (state == "PENDING") {
                out("R", spec, "pending", state, n, "", "waiting for a slot", "", "")
            } else if (state == "CANCELLED") {
                out("R", spec, "cancelled", state, n, "", "never started, cancelled" by, "", \
                    "CANCELLED before starting" by)
            } else {
                out("R", spec, "failed", state, n, "", "never started", "", state " before starting")
            }
            continue
        }

        t = id
        sub(/^[^_]*_/, "", t)
        see(t + 0)
        code = rexit[i]
        sub(/:.*/, "", code)
        if (state ~ /^(RUNNING|COMPLETING|CONFIGURING|SUSPENDED|STOPPED)$/) cat = "running"
        else if (state ~ /^(PENDING|REQUEUED|REQUEUE_HOLD|RESIZING)$/) cat = "pending"
        else if (state == "COMPLETED") cat = (banner[t] == "FAILED") ? "failed" : "succeeded"
        else if (state == "CANCELLED") cat = "cancelled"
        else cat = "failed"

        detail = ""
        pointer = ""
        reason = ""
        last = (lastproc[t] != "") ? lastproc[t] : lastline[t]
        if (!(t in have_log)) {
            detail = (cat == "running" || cat == "pending") ? "starting" : "(no log " job "_" t ".out)"
            if (cat == "failed") reason = state " · no log"
            if (cat == "cancelled") reason = "CANCELLED" (who != "" ? " by " who : "")
        } else if (cat == "succeeded") {
            detail = (counts[t] != "") ? counts[t] : "done"
        } else if (cat == "running") {
            if (retries[t] > 0) detail = retries[t] (retries[t] == 1 ? " retry" : " retries") " · "
            detail = detail "last: " (last != "" ? last : "starting")
        } else if (cat == "cancelled") {
            detail = (who != "" ? "by " who : "cancelled") (last != "" ? " · last: " last : "")
            pointer = base[t]
            reason = "CANCELLED" (who != "" ? " by " who : "")
        } else if (cat == "failed") {
            if (errproc[t] != "") {
                detail = errproc[t] " exit " errexit[t] exit_note(errexit[t])
                if (counts[t] != "") detail = detail " · " counts[t]
                reason = state " · nextflow " errproc[t] " exit " errexit[t]
                sub(/ \([0-9]+\)$/, "", reason)
            } else if (errmsg[t] != "") {
                detail = errmsg[t]
                reason = state " · " mask(errmsg[t])
            } else if (state == "TIMEOUT") {
                detail = "hit the walltime" (last != "" ? " · last: " last : "")
                reason = "TIMEOUT"
            } else if (state == "OUT_OF_MEMORY") {
                detail = "head process ran out of memory"
                reason = "OUT_OF_MEMORY (head process)"
            } else {
                detail = (lastline[t] != "") ? lastline[t] : "exit " rexit[i]
                reason = state " · " mask(detail)
            }
            pointer = base[t] (errline[t] != "" ? ":" errline[t] : "")
        }
        out("T", t, cat, state, relapsed[i], acc[t], detail, pointer, reason)
    }

    # pending tasks from squeue, as runs of consecutive tasks waiting for the same reason
    start = 0
    for (k = 1; k <= max_pend + 1; k++) {
        if (start > 0 && (k > max_pend || !(k in pend) || pend[k] != pend[start])) {
            out("R", range_spec(start, k - 1), "pending", "PENDING", k - start, "", \
                pending_detail(pend[start]), "", "")
            start = 0
        }
        if (start == 0 && k <= max_pend && (k in pend)) start = k
    }

    # A task split off the array to start next, and cancelled while still queued, leaves no
    # accounting record and no log. Report any index below the highest one seen that has no record
    # rather than silently dropping it from the totals.
    start = 0
    for (k = 1; k <= max_index + 1; k++) {
        if (k <= max_index && !(k in seen)) {
            if (start == 0) start = k
        } else if (start > 0) {
            out("R", range_spec(start, k - 1), "missing", "NO_RECORD", k - start, "", \
                "never started, and SLURM kept no record (likely queued when cancelled)", "", \
                "no SLURM record (never started)")
            start = 0
        }
    }
}
'

function collect_array {
    local task_logs=()
    while IFS= read -r task_log; do
        task_logs+=("$task_log")
    done < <(find "$log_folder" -maxdepth 1 -name "${job_id}_*.out" 2> /dev/null | sort -V)
    {
        sacct -j "$job_id" -X -n -P -o JobID,State,ExitCode,Elapsed | sed 's/^/S\t/'
        # Once the array has left the queue, squeue prints an error on stdout (see job_state), so
        # keep only task lines.
        squeue -r -h -j "$job_id" -o "%K|%T|%M|%r" 2> /dev/null \
            | { grep -E '^[0-9]+\|' || true; } \
            | sed 's/^/Q\t/' || true
        if [[ ${#task_logs[@]} -gt 0 ]]; then
            awk -v q="'" "$log_awk" ${task_logs[@]+"${task_logs[@]}"} | sed 's/^/L\t/'
        fi
    } | awk -v job="$job_id" "$join_awk"
}

# Draw the array summary from the rows collect_array produced.
# shellcheck disable=SC2016
render_awk='
function trunc(s) { return (length(s) > 110) ? substr(s, 1, 107) "..." : s }
function row(line,   f) {
    split(line, f, "\t")
    if (f[1] == "R") {
        printf "   %-27s %-13s %-10s %s\n", "tasks " f[2] " (" f[5] ")", f[4], "", f[7]
        return
    }
    printf "   %5s  %-20s %-13s %-10s %s\n", f[2], (f[6] == "" ? "-" : f[6]), f[4], f[5], trunc(f[7])
    if (f[8] != "") printf "   %5s  %-20s %-13s %-10s -> %s\n", "", "", "", "", f[8]
}
function header() {
    printf "   %5s  %-20s %-13s %-10s %s\n", "task", "accession", "state", "elapsed", "detail"
}
function section(title, list, n,   i, f, limit, tasks) {
    if (n == 0) return
    limit = (verbose == "true" || n <= max_rows) ? n : max_rows
    # a range row stands for several tasks, so count tasks rather than rows
    tasks = 0
    for (i = 1; i <= n; i++) {
        split(list[i], f, "\t")
        tasks += (f[1] == "R") ? f[5] : 1
    }
    printf "\n=== %s (%d)\n", title, tasks
    header()
    for (i = 1; i <= limit; i++) row(list[i])
    if (limit < n) printf "   ... and %d more rows (-v to list all)\n", n - limit
}
BEGIN { FS = "\t" }
{
    n = ($1 == "R") ? $5 + 0 : 1
    total += n
    cnt[$3] += n
    if ($3 == "failed" || $3 == "cancelled" || $3 == "missing") {
        # Keep failures ahead of cancellations, so a real error is not pushed past the row limit
        # (or down the reason list) by tasks that were cancelled on purpose.
        if ($3 == "failed") att_failed[++nf] = $0
        else if ($3 == "missing") att_missing[++nm] = $0
        else att_cancelled[++nc] = $0
        if ($9 != "") {
            if (!($9 in why)) {
                reasons[++nw] = $9
                rank[$9] = ($3 == "failed") ? 0 : ($3 == "missing") ? 1 : 2
            }
            why[$9] += n
        }
    } else if ($3 == "running") {
        run[++nr] = $0
    } else if ($3 == "succeeded") {
        ok[++ns] = $0
    } else if ($3 == "pending") {
        pend[++np] = $0
    }
}
END {
    printf "=== trickle array %s: %d tasks\n", job, total
    printf "    succeeded %d | failed %d | cancelled %d | running %d | pending %d%s\n", \
        cnt["succeeded"], cnt["failed"], cnt["cancelled"], cnt["running"], cnt["pending"], \
        (cnt["missing"] > 0 ? " | no record " cnt["missing"] : "")
    if (nw > 1) {
        # A handful of distinct reasons explains a long list, so show them first: failures, then
        # missing tasks, then cancellations, most common first within each (selection sort: there
        # are only ever a few).
        for (i = 1; i <= nw; i++)
            for (j = i + 1; j <= nw; j++)
                if (rank[reasons[j]] < rank[reasons[i]] \
                    || (rank[reasons[j]] == rank[reasons[i]] && why[reasons[j]] > why[reasons[i]])) {
                    tmp = reasons[i]
                    reasons[i] = reasons[j]
                    reasons[j] = tmp
                }
        printf "\n=== needs attention, by reason\n"
        for (i = 1; i <= nw; i++) printf "   %5d  %s\n", why[reasons[i]], reasons[i]
    }
    for (i = 1; i <= nf; i++) att[++na] = att_failed[i]
    for (i = 1; i <= nm; i++) att[++na] = att_missing[i]
    for (i = 1; i <= nc; i++) att[++na] = att_cancelled[i]
    section("needs attention", att, na)
    section("running", run, nr)
    if (verbose == "true") section("pending", pend, np)
    if (ns > 0) {
        if (verbose == "true") section("succeeded", ok, ns)
        else printf "\n=== succeeded: %d (-v to list)\n", ns
    }
    printf "\n   logs: %s/%s_*.out\n", folder, job
    printf "   one task in detail: %s %s_TASK\n", script, job
}
'

function render_array {
    local rows="$1"
    printf '%s\n' "$rows" \
        | awk -v job="$job_id" -v verbose="$verbose" -v max_rows="$max_rows" \
            -v folder="$log_folder" -v script="$0" "$render_awk"
    echo
    render_queue
}

function array_active {
    printf '%s\n' "$1" \
        | awk -F'\t' '$3 == "running" || $3 == "pending" {n += ($1 == "R") ? $5 : 1} END {print n + 0}'
}

function array_attention {
    # task, state, accession, detail for each task that failed or was cancelled
    printf '%s\n' "$1" \
        | awk -F'\t' '$1 == "T" && ($3 == "failed" || $3 == "cancelled") {print $2 "\t" $4 "\t" $6 "\t" $7}'
}

# A false positive in shellcheck 0.11.0 reports SC2218 ("only defined later") for each function
# called via $(...) in a function that also has another command substitution such as $(date).
# Reproduced in isolation; every function called here is defined above, so disable it here only.
# shellcheck disable=SC2218
function watch_array {
    local rows screen attention previous="" new events="" first=true stamp finished=false
    until [[ "$finished" == true ]]; do
        rows=$(collect_array)
        attention=$(array_attention "$rows")
        screen=$(render_array "$rows")
        new=""
        if [[ "$first" == false ]]; then
            # Tag each list so an empty previous list cannot be mistaken for the current one.
            new=$(
                {
                    printf '%s\n' "$previous" | sed 's/^/P\t/'
                    printf '%s\n' "$attention" | sed 's/^/C\t/'
                } | awk -F'\t' '$1 == "P" {seen[$2] = 1; next} $2 != "" && !($2 in seen)'
            )
            if [[ -n "$new" ]]; then
                stamp=$(date '+%H:%M:%S')
                events+=$(
                    printf '%s\n' "$new" \
                        | awk -F'\t' -v s="$stamp" \
                            '{printf "   %s  task %s  %s  %s: %s\n", s, $2, ($4 == "" ? "-" : $4), $3, $5}'
                )$'\n'
            fi
        fi
        clear_screen
        echo "status of array $job_id, refreshed every ${interval}s at $(date '+%F %T'); Ctrl-C to stop"
        if [[ -n "$events" ]]; then
            echo
            echo "=== newly failed or cancelled since watching started"
            printf '%s' "$events"
        fi
        echo
        printf '%s\n' "$screen"
        # only ring on a terminal: in a redirected log the bell is a stray byte
        if [[ -n "$new" && -t 1 ]]; then
            printf '\a'
        fi
        if [[ "$(array_active "$rows")" -eq 0 ]]; then
            echo
            echo "=== array $job_id has finished"
            finished=true
        else
            previous="$attention"
            first=false
            # The Sherlock login profile makes sleep a function that runs it under timeout, which
            # moves it out of the terminal's process group: Ctrl-C then never reaches it, and bash
            # waits for it regardless. command runs the real sleep, which Ctrl-C stops.
            command sleep "$interval"
        fi
    done
}

# A false positive in shellcheck 0.11.0 reports SC2218 ("only defined later") for each function
# called via $(...) in a function that also has another command substitution such as $(date).
# Reproduced in isolation; every function called here is defined above, so disable it here only.
# shellcheck disable=SC2218
function watch_single {
    local screen state finished=false
    until [[ "$finished" == true ]]; do
        screen=$(render_single)
        state=$(job_state "$job_id")
        clear_screen
        echo "status of job $job_id, refreshed every ${interval}s at $(date '+%F %T'); Ctrl-C to stop"
        echo
        printf '%s\n' "$screen"
        if [[ -z "$state" ]]; then
            echo
            echo "=== job $job_id has left the queue"
            finished=true
        else
            # command: see watch_array
            command sleep "$interval"
        fi
    done
}

if [[ "$is_array" == true ]]; then
    if [[ "$watch" == true ]]; then
        watch_array
        exit 0
    fi
    render_array "$(collect_array)"
    if [[ "$follow" == true ]]; then
        task_logs=()
        while IFS= read -r task_log; do
            task_logs+=("$task_log")
        done < <(find "$log_folder" -maxdepth 1 -name "${job_id}_*.out" 2> /dev/null | sort -V)
        if [[ ${#task_logs[@]} -gt 0 ]]; then
            echo
            echo "=== following ${#task_logs[@]} task logs under $log_folder"
            tail -f ${task_logs[@]+"${task_logs[@]}"}
        fi
    fi
    exit 0
fi

if [[ "$watch" == true ]]; then
    watch_single
    exit 0
fi
render_single
if [[ "$follow" == true && -f "$log_folder/$job_id.out" ]]; then
    echo
    echo "=== following $log_folder/$job_id.out"
    tail -f "$log_folder/$job_id.out"
fi
