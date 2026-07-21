#!/usr/bin/env bash
set -euo pipefail
# script to build a docker image for a conda yaml, using pixi as the builder and maintaining the
# lock file. Note that this would not work for an actual pixi project because no source is copied,
# only the .toml and .lock

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
yaml=$1
tag="${2:-latest}"
registry="${3:-kundajelab}"

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
repo_dir=$(dirname "$script_dir")
dockerfile="$repo_dir/dockers/pixi-yaml.Dockerfile"
pushd "$repo_dir" &> /dev/null

# first get .lock and .toml, work in temporary folder inside environment folder, so they're on the
# same device and we can use hardlinks (pixi won't follow symlinks)
environment_name=$(basename "${yaml%.yaml}")
lock="${yaml%.yaml}.lock"
toml="${yaml%.yaml}.toml"
temp_dir="$(mktemp -d -p "$(dirname "$yaml")")"
trap 'rm -rf "$temp_dir"' EXIT
cp "$yaml" "$temp_dir/"
if [[ -f "$lock" ]]; then
    # lock already exists. Use existing lock for faster update
    ln "$lock" "$temp_dir/pixi.lock"
fi
# create .toml from scratch
pushd "$temp_dir" &> /dev/null
unset PIXI_PROJECT_MANIFEST
pixi init \
    --import "${environment_name}.yaml" \
    -p linux-64 -p osx-64 -p osx-arm64 -p linux-aarch64 \
    .
# update lock
pixi lock --no-install --manifest-path "./pixi.toml"
popd &> /dev/null
cp "$temp_dir/pixi.toml" "$toml"
if [[ ! -f "$lock" ]]; then
    cp "$temp_dir/pixi.lock" "$lock"
fi

# now build docker
pushd "$(dirname "$yaml")" &> /dev/null
container_name=$(tr '[:upper:]' '[:lower:]' <<< "$environment_name")
local_tag="${container_name}:${tag}"
remote_tag="${registry}/${local_tag}"
docker build \
    --platform linux/amd64,linux/arm64 \
    --build-arg ENV_NAME="$environment_name" \
    -t "${local_tag}" \
    -t "${remote_tag}" \
    -f "$dockerfile" \
    .

popd &> /dev/null

if [[ "$push" == "true" ]]; then
    # push to remote repo
    1>&2 docker push "${remote_tag}"
fi

# output name for dotenv environment file
echo "${environment_name}_IMAGE=${remote_tag}"
