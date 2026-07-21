#!/usr/bin/env bash
set -euo pipefail

args=()
while [[ "$#" -ge 1 ]]; do
    case "$1" in
        "--"*)
            args+=("$1")
            shift 1
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

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
repo_dir=$(dirname "$script_dir")

"$script_dir/find-projects.sh" \
| while read -r project type; do
    case "$type" in
        "pixi")
            1>&2 echo "Building docker for pixi project: $project"
            "$script_dir/build-project-docker.sh" "${args[@]}" "$project" "$@"
            ;;
        "uv")
            1>&2 echo "Building docker for uv project: $project"
            "$script_dir/build-project-docker.sh" "${args[@]}" "$project" "$@"
            ;;
        "yaml")
            1>&2 echo "Building docker for yaml project: $project"
            "$script_dir/build-yaml-docker.sh" "${args[@]}" "$repo_dir/environments/$project" "$@"
            ;;
        *)
            1>&2 echo "Unknown project type: $type for project: $project"
            exit 1
            ;;
    esac
done \
> "$repo_dir/.env"