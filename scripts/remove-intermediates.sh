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

ARGS:
    -h|--help: show this message and exit.
    -w|--workspace: Use this as the root workspace instead of default.

EOF
}

workspace="$("$script_dir/get-default-workspace.sh")"
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
        "--")
            break
            shift 1
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
    if [[ "$(basename "$project_dir")" == "apptainer_cache" ]]; then
        return  # don't clear the common apptainer_cache
    elif [[ ! -d "$project_dir/output" ]]; then
        # this is not a project folder, remove the whole thing
        rm -r "$project_dir"
        return
    fi

    find "$project_dir/output" -type l \
    | while read -r link_path; do
        fix_link "$link_path"
    done

    rm -r "$project_dir/work"
}

if [[ "$1" == "all" ]]; then
    find "$workspace" -mindepth 1 -maxdepth 1 -type d \
        | while read -r project_dir; do
            if [[ "$(basename "$project_dir")" == "apptainer_cache" ]]; then
                continue  # don't clear the common apptainer_cache
            fi
            clear_workspace "$project_dir"
        done
    rm -rf "$repo_dir/.nextflow*"
else
    clear_workspace "$workspace/$1"
fi
