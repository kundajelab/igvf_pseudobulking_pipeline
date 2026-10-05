#!/usr/bin/env bash
set -euo pipefail

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
repo_dir=$(dirname "$script_dir")

function usage {
    cat << EOF
Usage: $0 [ARGS] -- ANALYSIS_ACCESSION

Remove all intermediate files in the project output folder, replacing symlinks in the output with
regular files.

If ANALYSIS_ACCESSION is "all", do this for every project folder.

A project folder whose .nextflow.log changed in the last ACTIVE_MINUTES is assumed to have a run in
progress, and is skipped.

ARGS:
    -h|--help: show this message and exit.
    -w|--workspace: Use this as the root workspace instead of default.
    -a|--active-minutes: Skip folders whose nextflow log changed this recently. Default: $active_minutes

EOF
}

workspace="$("$script_dir/get-default-workspace.sh")"
# NOTE: while tasks are running, nextflow writes its pending-task list to the log every 5 minutes
active_minutes=30
while [[ "$#" -ge 1 ]]; do
    case "$1" in
        "-h" | "--help")
            usage
            exit 0
            ;;
        "-w" | "--workspace")
            workspace="$2"
            shift 2
            ;;
        "-a" | "--active-minutes")
            active_minutes="$2"
            shift 2
            ;;
        "--")
            shift 1
            break
            ;;
        *)
            break
            ;;
    esac
done

function fix_link {
    local -r link_path="$1"
    original_path=$(readlink -f "$link_path")
    rm "$link_path"
    cp -f "$original_path" "$link_path"
}

function clear_workspace {
    local -r project_dir="$1"
    case "$(basename "$project_dir")" in
        apptainer_cache|conda_cache)
            return  # don't clear the common container and environment caches
            ;;
    esac
    # run-pipeline.sh launches nextflow from inside the project folder, so its log is written here
    if [[ -n "$(find "$project_dir" -maxdepth 1 -name .nextflow.log -mmin "-$active_minutes")" ]]; then
        1>&2 echo "Skipping $project_dir: its nextflow log changed in the last $active_minutes minutes, so a run may be in progress"
        return
    fi
    1>&2 echo "Clearing $project_dir"
    if [[ ! -d "$project_dir/output" ]]; then
        # this is not a project folder, remove the whole thing
        rm -rf "$project_dir"
        return
    fi

    find "$project_dir/output" -type l \
    | while read -r link_path; do
        fix_link "$link_path"
    done

    rm -rf "$project_dir/work"
}

if [[ "$#" -lt 1 ]]; then
    1>&2 usage
    exit 1
fi

if [[ "$1" == "all" ]]; then
    find "$workspace" -mindepth 1 -maxdepth 1 -type d \
        | while read -r project_dir; do
            clear_workspace "$project_dir"
        done
    rm -rf "$repo_dir/.nextflow"*
else
    clear_workspace "$workspace/$1"
fi
