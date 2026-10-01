#!/usr/bin/env bash
set -euo pipefail

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

#tag="$(git describe --tags --abbrev=0)"
project="all"
workspace="$(scripts/get-default-workspace.sh)"
apptainer_cache_dir="$workspace/apptainer_cache"

function usage {
    cat << EOF
Usage: $0 [ARGS] -- command [command args]

Run the specified command in the specified project container image.
At the moment, only works with apptainer, not docker.
If no command is specified, just pull the image.
To run interactively pass the bash command (or some other shell).

ARGS:
    -h|--help: Show this message and exit.
    -p|--project: project name for apptainer image
    -t|--tag: docker/apptainer image tag. Defaults to the tag in the repo .env file.
    -c|--cache: Where apptainer images are cached. Defaults to inferred from environment ($apptainer_cache_dir)
EOF
}

while [[ "$#" -ge 1 ]]; do
    case "$1" in
        "-h" | "--help")
            usage
            exit 0
            ;;
        "-p" | "--project")
            project="$2"
            shift 2
            ;;
        "-t" | "--tag")
            tag="$2"
            shift 2
            ;;
        "-c" | "--cache")
            apptainer_cache_dir="$2"
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

function pull {
    # first pull the image with minimal nuisance warnings
    local -r project="$1"
    1>&2 echo "Pulling $project"
    if [[ -z "${tag+x}" ]]; then
        # tag isn't defined, get it from .env
        local -r pull_tag="$("$script_dir/get-project-tag.sh" "$project")"
    else
        local -r pull_tag="$tag"
    fi

    if [[ "$pull_tag" =~ sha256 ]]; then
        local -r apptainer_img="$apptainer_cache_dir/kundajelab-${project}@${pull_tag//:/-}.img"
        local -r remote_img="docker://kundajelab/${project}@$pull_tag"
    else
        local -r apptainer_img="$apptainer_cache_dir/kundajelab-${project}-$pull_tag.img"
        local -r remote_img="docker://kundajelab/${project}:$pull_tag"
    fi
    apptainer pull \
        -F \
        "$apptainer_img" \
        "$remote_img" \
        1>&2 \
        2> >(grep -v "harmless EPERM" >&2)

    echo "$apptainer_img"
}

if [[ "$project" == "all" ]]; then
    # project set to all
    if [[ "$#" -gt 0 ]]; then
        1>&2 echo "Must specify a project to run a command."
        1>&2 usage
        exit 1
    fi
    # no command, just pull / update all projects
    script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
    while read -r project _project_type; do
        pull "$project"
    done < <("$script_dir/find-projects.sh")
else
    apptainer_img=$(pull "$project")

    if [[ "$#" -gt 0 ]]; then
        command="${*}"
        apptainer run "$apptainer_img" "${command[@]}"
    fi
fi
