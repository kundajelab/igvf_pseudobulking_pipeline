#!/usr/bin/env bash
set -euo pipefail

tag="$(git describe --tags --abbrev=0)"
project=""
workspace="$(scripts/get-default-workspace.sh)"
apptainer_cache_dir="$workspace/apptainer_cache/"

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
    -t|--tag: docker/apptainer image tag ($tag)
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

# first pull the image with minimal nuisance warnings
apptainer pull \
    -F \
    "$apptainer_cache_dir/kundajelab-${project}-$tag.img" \
    "docker://kundajelab/${project}:$tag" \
    2> >(grep -v "harmless EPERM" >&2)

if [[ "$#" -gt 0 ]]; then
    command="${*}"
    apptainer run "$apptainer_cache_dir/kundajelab-${project}-$tag.img" "${command[@]}"
fi
