#!/usr/bin/env bash
set -euo pipefail

function usage {
    cat << EOF
Usage: $0 [ARGS] -- update_command [update_value]

Update the project version, as well as sub-projects. Tag git with version.

If the update_command is:
- set: then the update_value must be specified and a valid semantic version string
  ([major].[minor].[patch] see https://semver.org).
- major: bump the major version number.
- minor: bump the minor version number.
- patch: bump the patch version number.

ARGS:
    -h|--help: show this message and exit.
    -m|--message: Set message for git tag.
    -c|--commit: Which commit to tag (defaults to HEAD)
    --push/--no-push: Whether or not to push the git tag to origin

EOF
}

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
repo_dir=$(dirname "$script_dir")
pushd "$repo_dir" &> /dev/null

commit=HEAD
# determine if we are setting or bumping the project version
push="false"
while [[ "$#" -ge 1 ]]; do
    case "$1" in
        "-h" | "--help")
            usage
            exit 0
            ;;
        "-m" | "--message")
            message="$2"
            shift 2
            ;;
        "-c" | "--commit")
            commit="$2"
            shift 2
            ;;
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

function validate_version {
    version="$1"
    # validate the version string
    if [[ ! "$version" =~ ^v?[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
        1>&2 echo "Invalid version: '$version'"
        1>&2 echo "Version must adhere to semantic versioning: https://semver.org"
        1>&2 echo "[major].[minor].[patch]" with major, minor, and patch being integers.
        exit 1
    fi
    # remove leading v if present
    sed 's/^v//' <<< "$version"
}

command="$1"
case "$command" in
    "set")
        version="$(validate_version "$2")"
        pixi workspace version set --manifest-path . "$version"
        ;;
    "major"|"minor"|"patch")
        pixi workspace version "$command"
        version="$(pixi workspace version get --manifest-path .)"
        ;;
    *)
        1>&2 echo "Invalid version adjustment command: '$command'"
        exit 1
        ;;
esac

if [[ -z "${message:-}" ]]; then
    message="Released version $version"
fi

"$script_dir/find-projects.sh" \
| while read -r project type; do
    case "$type" in
        "pixi")
            1>&2 echo "Updating version for pixi project: $project"
            pixi workspace version set --manifest-path "$project" "$version"
            ;;
        "uv")
            1>&2 echo "Updating version for uv project: $project"
            pushd "$repo_dir/$project" &> /dev/null
            pixi run uv version "$version"
            popd &> /dev/null
            ;;
        "yaml")
            # Nothing to do
            ;;
        *)
            1>&2 echo "Unknown project type: $type for project: $project"
            exit 1
            ;;
    esac
done

git commit -m "${message}"
git tag \
    -a "v${version}" \
    -m "${message}" \
    "${commit}"
if [[ "$push" == "true" ]]; then
    git push origin "v${version}"
fi
