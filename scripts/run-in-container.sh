#!/usr/bin/env bash
set -euo pipefail

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

project="all"
workspace="$("$script_dir/get-default-workspace.sh")"
apptainer_cache_dir="$workspace/apptainer_cache"
work_dir=$(pwd)

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
    -q|--queue: If specified, submit the command to run in that queue. If empty, execute normally. (empty)
    -W|--work-dir: Where to execute the command. Defaults to pwd ($work_dir)
    -T|--time|-C|--cpus|-M|--mem|-l|--log-dir|-n|--name": Pass to submit-command.sh if queue is specified.
EOF
}

queue=""
submit_args=()
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
        "-q" | "--queue")
            queue="$2"
            shift 2
            ;;
        "-W" | "--work-dir")
            work_dir="$2"
            shift 2
            ;;
        "-T" | "--time" | "-C" | "--cpus" | "-M" | "--mem" | "-l" | "--log-dir" | "-n" | "--name")
            submit_args+=("$1" "$2")
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
    while read -r project _project_type; do
        pull "$project"
    done < <("$script_dir/find-projects.sh")
else
    apptainer_img=$(pull "$project")

    if [[ "$#" -gt 0 ]]; then
        if [[ -n "$queue" ]]; then
            "$script_dir/submit-command.sh" \
                --queue "$queue" \
                --name "$(IFS=- ; echo "$*")" \
                --work-dir "$work_dir" \
                ${submit_args[@]+"${submit_args[@]}"} \
                -- \
                apptainer run "$apptainer_img" "$@"
        elif [[ "${#submit_args[@]}" -gt 0 ]]; then
            1>&2 echo "submit-args were specified without queue"
            usage
            exit 1
        else
            pushd "$work_dir" &> /dev/null
            apptainer run "$apptainer_img" "$@"
        fi
    fi
fi
