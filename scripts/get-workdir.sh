#!/usr/bin/env bash
set -euo pipefail

workspace="$1"
job_id="$2"
work_dir=$(find "${workspace}/work/$(dirname "$job_id")" -maxdepth 1 -type d -name "$(basename "$job_id")*" | tail -n1)
1>&2 echo "work_dir: ${work_dir}"
echo "${work_dir}"
