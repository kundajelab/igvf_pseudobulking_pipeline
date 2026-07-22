#!/usr/bin/env bash
set -euo pipefail

workspace="$1"
job_id="$2"
first_dir=$(dirname "$job_id")
second_dir=$(basename "$job_id")
work_dir=$(\
    find workspace -type d -mindepth 2 -maxdepth 2 -name work\
    | while read -r workdir; do
        sub_dir="$workdir/$first_dir"
        if [[ -d "$sub_dir" ]]; then
            find "$sub_dir" -maxdepth 1 -type d -name "${second_dir}*"
        fi
    done \
    | tail -n1
)
1>&2 echo "work_dir: ${work_dir}"
echo "${work_dir}"
