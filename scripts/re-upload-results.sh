#!/usr/bin/env bash
set -euxo pipefail

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
repo_dir=$(dirname "$script_dir")

queue="akundaje"
workspace="$("$script_dir/get-default-workspace.sh")"
mode="prod"
dry_run_arg="--no-dry-run"

function usage {
    cat << EOF
Usage: $(basename "$0") [ARGS] metadata

Re-upload completed pipeline results in a SLURM batch job.

metadata is specified by accession ID.

ARGS:
    -h|--help: Show this message and exit.
    -q|--queue: SLURM, use this slurm queue. Default: $queue
    -w|--workspace: Where to output files. Defaults to inferred from environment ($workspace)
    -T|--time|-C|--cpus|-M|--mem|-l|--log-dir|-n|--name": Pass to submit-command.sh if queue is specified.
    -m|--mode: Specify IGVF server: "prod" or "staging" ($mode)
    --dry-run/--no-dry-run: Turn on/off igvf_dry_run. With dry-run *on* no changes to the IGVF portal are made.
      With dry-run *off* records are posted/patched and files are uploaded. Default: $dry_run_arg
EOF
}
submit_args=()

while [[ "$#" -ge 1 ]]; do
    case "$1" in
        "-h" | "--help")
            usage
            exit 0
            ;;
        "-T" | "--time" | "-C" | "--cpus" | "-M" | "--mem" | "-l" | "--log-dir" | "-n" | "--name")
            submit_args+=("$1" "$2")
            shift 2
            ;;
        "-w" | "--workspace")
            workspace="$2"
            shift 2
            ;;
        "-q" | "--queue")
            queue="$2"
            shift 2
            ;;
        "-m" | "--mode")
            mode="$2"
            case "$mode" in
                prod|staging)
                    ;;
                *)
                    1>&2 echo "Invalid IGVF portal mode: $mode"
                    exit 1
                    ;;
            esac
            shift 2
            ;;
        "--dry-run"|"--no-dry-run")
            dry_run_arg="$1"
            shift 1
            ;;
        --)
            shift 1
            break
            ;;
        "--"?*)
            1>&2 echo "Unknown argument: $1"
            exit 1
            ;;
        *)
            break
            ;;
    esac
done

if [[ $# != 1 ]]; then
    1>&2 echo "Metadata accession must be specified."
    usage
    exit 1
fi
if ! grep -q IGVF_API_KEY <<< "$(env)" || ! grep -q IGVF_SECRET_KEY <<< "$(env)"; then
    1>&2 echo "IGVF_API_KEY and IGVF_SECRET_KEY must be set."
    exit 1
fi

metadata_accession="$1"
metadata_file="${metadata_accession}.tsv.gz"

# Find the run folder from the metadata that run-pipeline.sh downloaded into it, rather than
# asking the portal, which would need to know which portal the dry run was made for.
mapfile -t run_folders < <(
    find "$workspace" -mindepth 2 -maxdepth 2 -type f -name "$metadata_file" -exec dirname {} \;
)
if [[ "${#run_folders[@]}" -ne 1 ]]; then
    1>&2 echo "Expected one run folder in '$workspace' holding '$metadata_file', found ${#run_folders[@]}."
    exit 1
fi
run_folder="${run_folders[0]}"
analysis_accession="$(basename "$run_folder")"

# find the working folder for the IGVF_UPLOAD process
pushd "$run_folder" &> /dev/null
# NOTE: output/upload.sh is published as a symlink into the IGVF_UPLOAD work dir. If it is missing,
# or is a regular file because remove-intermediates.sh replaced the symlinks, readlink -f still
# succeeds and would point at the output folder, which lacks the staged inputs.
if [[ ! -L output/upload.sh ]]; then
    1>&2 echo "'$run_folder/output/upload.sh' is not a symlink into the IGVF_UPLOAD work dir."
    exit 1
fi
# NOTE: resolve the link on its own line. For a dangling link, readlink -f prints nothing and fails,
# which a nested $(dirname "$(readlink -f ...)") would hide, yielding "." (the run folder).
upload_script="$(readlink -f output/upload.sh || true)"
if [[ ! -f "$upload_script" ]]; then
    1>&2 echo "'$run_folder/output/upload.sh' links to a missing file; was the work dir removed?"
    exit 1
fi
scratch_dir="$(dirname "$upload_script")"
if [[ ! -e "$scratch_dir/$metadata_file" ]]; then
    1>&2 echo "IGVF_UPLOAD work dir '$scratch_dir' does not hold '$metadata_file'."
    exit 1
fi
pushd "$scratch_dir" &> /dev/null

# write a fix script there
cat << EOF > fix.sh
#!/usr/bin/env bash
set -euo pipefail

# re-write the upload script
pushd "$scratch_dir" &> /dev/null

igvf-portal gen-upload-script \
    . \
    --input-file-sets "$analysis_accession" \
    --metadata-file "$metadata_file" \
    --annotations-tsv "./cell_name_to_annotation_mapping.tsv"\
    "$dry_run_arg" \
    --igvf-mode "$mode"

./upload.sh
EOF

chmod u+x fix.sh

# run the fix script in the igvf_portal image on the requested queue
pixi run --manifest-path "$repo_dir" \
    run-in-container \
    --project igvf_portal \
    ${submit_args[@]+"${submit_args[@]}"} \
    --work-dir "$scratch_dir" \
    --queue "$queue" \
    -- \
    ./fix.sh
