#!/usr/bin/env bash
set -euo pipefail

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
repo_dir=$(dirname "$script_dir")

"$script_dir/find-projects.sh" \
| while read -r project type; do
    case "$type" in
        "pixi")
            echo "Installing pixi project: $project"
            pushd "$repo_dir/$project" &> /dev/null
            pixi install --manifest-path .
            popd &> /dev/null
            ;;
        "uv")
            echo "Installing uv project: $project"
            pushd "$repo_dir/$project" &> /dev/null
            pixi run uv build && pixi run uv sync
            popd &> /dev/null
            ;;
        "yaml")
            # no need to install yaml projects
            ;;
        *)
            echo "Unknown project type: $type for project: $project"
            exit 1
            ;;
    esac
done