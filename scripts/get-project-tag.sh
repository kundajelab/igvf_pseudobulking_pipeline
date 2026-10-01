#!/usr/bin/env bash
set -euo pipefail

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
repo_dir=$(dirname "$script_dir")

project="$1"
environment_name=$(tr '[:lower:]' '[:upper:]' <<< "$project")
env_pattern="${environment_name}_IMAGE=kundajelab/${project}"
grep "^${env_pattern}" "$repo_dir/.env" \
    | cut -c "$((2+${#env_pattern}))-"
