#!/usr/bin/env bash
set -euo pipefail

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
repo_dir=$(dirname "$script_dir")

find "$repo_dir" -mindepth 1 -maxdepth 1 -type d \
| while read -r dir; do
    if [[ -d "$dir/.pixi" ]]; then
        printf "%s\tpixi\n" "$(basename "$dir")"
    elif [[ -d "$dir/.venv" ]]; then
        printf "%s\tuv\n" "$(basename "$dir")"
    fi
done

find "$repo_dir/environments" -maxdepth 1 -type f -name "*.yaml" \
| while read -r yaml; do
    if [[ ! -d "$repo_dir/$(basename "${yaml%.yaml}" | tr '[:upper:]' '[:lower:]')" ]]; then
        printf "%s\tyaml\n" "$(basename "$yaml")"
    fi
done
