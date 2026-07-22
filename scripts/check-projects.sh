#!/usr/bin/env bash
set -euo pipefail

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
repo_dir=$(dirname "$script_dir")

"$script_dir/find-projects.sh" \
| while read -r project type; do
    case "$type" in
        "pixi")
            echo "Checking pixi project: $project"
            pushd "$repo_dir/$project" &> /dev/null
            pixi run --manifest-path . checks
            popd &> /dev/null
            ;;
        "uv")
            echo "Checking uv project: $project"
            pushd "$repo_dir/$project" &> /dev/null
            pixi run poe checks
            popd &> /dev/null
            ;;
        "yaml")
            echo "Checking yaml project: $project"
            environments_dir="$repo_dir/environments"
            yaml="$environments_dir/$(tr '[:upper:]' '[:lower:]' <<< "$project").yaml"
            yamllint \
                -d "{extends: default, rules: {document-start: false}}" \
                "$yaml"
            has_needed_keys=$(yq 'keys | contains(["channels", "dependencies", "name"])' "$yaml")
            if [[ "$has_needed_keys" != "true" ]]; then
                1>&2 echo "Error: $yaml is missing one or more required keys: channels, dependencies, name"
                exit 1
            fi
            ;;
        *)
            1>&2 echo "Unknown project type: $type for project: $project"
            exit 1
            ;;
    esac
done
