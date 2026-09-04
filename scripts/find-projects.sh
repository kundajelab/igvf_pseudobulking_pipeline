#!/usr/bin/env bash
set -euo pipefail

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
repo_dir=$(dirname "$script_dir")

find "$repo_dir" -mindepth 1 -maxdepth 1 -type d \
| while read -r dir; do
    if [[ -f "$dir/pixi.lock" ]]; then
        printf "%s\tpixi\n" "$(basename "$dir")"
    elif [[ -f "$dir/uv.lock" ]]; then
        printf "%s\tuv\n" "$(basename "$dir")"
    fi
done

find "$repo_dir/environments" -maxdepth 1 -type f -name "*.yaml" \
| while read -r yaml; do
    project=$(basename "${yaml%.yaml}" | tr '[:upper:]' '[:lower:]')
    if [[ ! -d "$repo_dir/$project" ]]; then
        printf "%s\tyaml\n" "$project"
    fi
done
