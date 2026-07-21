#!/usr/bin/env bash
set -euo pipefail
# script to build a docker image for a uv or pixi-based project

push="true"
while [[ "$#" -gt 1 ]]; do
    case "$1" in
        "--push")
            push="true"
            shift 1
            ;;
        "--no-push")
            push="false"
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
project="$1"
tag="${2:-latest}"
registry="${3:-kundajelab}"

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
repo_dir=$(dirname "$script_dir")
pushd "$project" &> /dev/null
if [[ -f "$repo_dir/$project/pixi.lock" ]]; then
    project_type="pixi"
    pixi lock --manifest-path .
else
    project_type="uv"
    pixi run uv build
    pixi run uv sync
fi
dockerfile="$repo_dir/dockers/${project_type}-project.Dockerfile"


local_tag="${project}:${tag}"
remote_tag="${registry}/${local_tag}"

docker build \
    -D \
    --platform linux/amd64,linux/arm64 \
    --build-arg ENV_NAME="$project" \
    -t "${local_tag}" \
    -t "${remote_tag}" \
    -f "$dockerfile" \
    .

# push to remote repo
if [[ "$push" == "true" ]]; then
    1>&2 docker push "${remote_tag}"
fi

# output name for dotenv environment file
environment_name=$(tr '[:lower:]' '[:upper:]' <<< "$project")
echo "${environment_name}_IMAGE=${remote_tag}"

