#!/usr/bin/env bash
set -euo pipefail

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
repo_dir=$(dirname "$script_dir")
pushd &> /dev/null "$repo_dir"

# The test metadata and the accessions it refers to, matching submit-pipeline.sh. Nothing is
# downloaded or submitted, but the required params have to be present for the config to parse and for
# validateParameters() to pass.
metadata_file="$repo_dir/test_metadata.tsv"
principal_analysis="IGVFDS5417HJRJ,IGVFDS6430MYNQ"
profile="$(scripts/get-default-profile.sh)"

function usage {
    cat << EOF
Usage: $0 [ARGS]

Build the whole workflow with "nextflow run -preview": every channel is wired up and every process
directive is resolved, but no task is submitted. This catches wiring mistakes that "nextflow lint"
cannot see, because the linter checks syntax without resolving channel operators or process inputs.

Caught by this, and NOT by nextflow lint:
  * a channel operator applied to something that is not a channel, e.g. ".buffer()" on a plain List
  * a per-task directive that reads an input without being a closure, e.g. "queue dry_run ? a : b"

NOT caught by either, and still needing a real run:
  * a misspelled channel operator such as ".ifempty" instead of ".ifEmpty". Nextflow only resolves
    the name when the upstream channel terminates, which never happens while no task runs.

ARGS:
    -h|--help: Show this message and exit.
    -p|--profile: Comma-separated nextflow profiles. Default: $profile
    -f|--metadata-file: Metadata TSV to preview with. Default: $metadata_file
    -a|--principal-analysis: Principal analysis accessions. Default: $principal_analysis
EOF
}

while [[ "$#" -ge 1 ]]; do
    case "$1" in
        "-h" | "--help")
            usage
            exit 0
            ;;
        "-p" | "--profile")
            profile="$2"
            shift 2
            ;;
        "-f" | "--metadata-file")
            metadata_file="$2"
            shift 2
            ;;
        "-a" | "--principal-analysis")
            principal_analysis="$2"
            shift 2
            ;;
        *)
            1>&2 echo "Unknown argument: $1"
            exit 1
            ;;
    esac
done

# Send the work directory, the trace files and the nextflow log to a scratch directory, so a preview
# never disturbs a real run's workspace or resume cache.
# NOTE: "-log" is a nextflow option rather than a "run" option, so it has to come before "run".
temp_dir=$(mktemp -d)
trap 'rm -rf "$temp_dir"' EXIT

# NOTE: species is deliberately left unset so that the preview exercises the DETERMINE_SPECIES branch
# of the workflow, which is the more complicated of the two. No portal lookup happens, because
# -preview does not run the process.
NXF_ANSI_LOG=false nextflow -log "$temp_dir/.nextflow.log" run . \
    -preview \
    -profile "$profile" \
    --metadata_file "$metadata_file" \
    --principal_analysis "$principal_analysis" \
    --workspace "$temp_dir/workspace"
