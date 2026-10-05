#!/usr/bin/env bash
set -euo pipefail

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
repo_dir=$(dirname "$script_dir")
pushd &> /dev/null "$repo_dir"

queue="normal"
profile="$(scripts/get-default-profile.sh)"
workspace="$(scripts/get-default-workspace.sh)"
time_limit="1-00:00:00"
cpus=2
mem="8G"
log_dir="$HOME/logs"

function usage {
    cat << EOF
Usage: $0 [ARGS] -- [metadata]

For a pipeline that succeeded with --dry-run for uploads to the IGVF portal, actually perform the upload.
Note: this uses scripts/run-in-container.sh which currently only works with apptainer, so it must be run
on sherlock.

The upload goes to the IGVF portal (prod or staging) that the dry run was generated for, because the files
it refers to only exist on that one.

ARGS:
    -h|--help: Show this message and exit.
    -p|--profile: Use comma-separated nextflow profiles. Defaults to inferred from environment ($profile)
    -q|--queue: If running via SLURM, use this queue. Default: $queue
    -w|--workspace: Where to output files. Defaults to inferred from environment ($workspace)
    -a|--principal-analysis: Specify the accession of the principal analysis set. Otherwise it is found
      from the run folder in the workspace holding the downloaded metadata.
EOF
}

principal_analysis=""
while [[ "$#" -ge 1 ]]; do
    case "$1" in
        "-h" | "--help")
            usage
            exit 0
            ;;
        "-T" | "--time")
            time_limit="$2"
            shift 2
            ;;
        "-C" | "--cpus")
            cpus="$2"
            shift 2
            ;;
        "-M" | "--mem")
            mem="$2"
            shift 2
            ;;
        "-p" | "--profile")
            profile="$2"
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
        "-a" | "--principal-analysis")
            principal_analysis="$2"
            shift 2
            ;;
        "--")
            shift 1
            break
            ;;
        "-"?*)
            1>&2 echo "Unknown argument: $1"
            exit 1
            ;;
        *)
            break
            ;;
    esac
done

pushd "$repo_dir" &> /dev/null

metadata="$1"
if [[ "$metadata" =~ \.tsv(\.gz)?$ ]]; then
    metadata_file="$metadata"
    # ensure we have the principal analysis accession
    if [[ -z "$principal_analysis" ]]; then
        1>&2 echo "Must specify metadata accession, or metadata file and principal analysis accession"
        exit 1
    fi
    run_folder="$workspace/${principal_analysis//,/-}"
else
    if [[ -n "$principal_analysis" ]]; then
        run_folder="$workspace/${principal_analysis//,/-}"
    else
        # Find the run folder from the metadata that run-pipeline.sh downloaded into it, rather than
        # asking the portal, which would need to know which portal the dry run was made for.
        mapfile -t run_folders < <(
            find "$workspace" -mindepth 2 -maxdepth 2 -type f -name "${metadata}.tsv.gz" -exec dirname {} \;
        )
        if [[ "${#run_folders[@]}" -ne 1 ]]; then
            1>&2 echo "Expected one run folder in '$workspace' holding '${metadata}.tsv.gz', found ${#run_folders[@]}."
            1>&2 echo "Specify the principal analysis with --principal-analysis."
            exit 1
        fi
        run_folder="${run_folders[0]}"
        # NOTE: the folder name has "-" in place of the "," between multiple accessions
        principal_analysis=$(basename "$run_folder")
    fi
    metadata_file="$run_folder/${metadata}.tsv.gz"
fi
if [[ ! -f "$metadata_file" ]]; then
    1>&2 echo "Unable to find metadata file '$metadata_file'."
    exit 1
fi

output_folder="$run_folder/output"
dry_run_upload_script="$output_folder/upload.sh"
earnest_upload_script="$output_folder/upload-no-dry-run.sh"
if [[ ! -f "$dry_run_upload_script" ]]; then
    1>&2 echo "Unable to find dry-run upload script '$dry_run_upload_script'."
    exit 1
fi
# The upload script targets the portal it was generated for, and the files it uploads only exist on
# that one. Report which it is.
mode=$(sed -n 's/^igvf_mode="\(.*\)"$/\1/p' "$dry_run_upload_script")
if [[ -z "$mode" ]]; then
    1>&2 echo "Unable to find the IGVF portal mode in '$dry_run_upload_script'."
    exit 1
fi
sed 's/^dry_run_arg=".*"$/dry_run_arg=""/' "$dry_run_upload_script" > "$earnest_upload_script"
chmod u+x "$earnest_upload_script"

# Need to run upload in a non-preemptible queue, not the login node
job_name=igvf-upload/$principal_analysis
log_folder="$log_dir/$job_name"
mkdir -p "$log_folder"
sbatch_script=$(cat << EOF
#!/usr/bin/env bash
#SBATCH --job-name=$job_name
#SBATCH --output=$log_folder/%j.out
#SBATCH --partition=$queue
#SBATCH --time=$time_limit
#SBATCH --cpus-per-task=$cpus
#SBATCH --mem=$mem
#SBATCH --chdir=$repo_dir
set -euo pipefail

echo "head process: job \$SLURM_JOB_ID on \$(hostname), partition \$SLURM_JOB_PARTITION"

pixi run run-in-container --project igvf_portal "$earnest_upload_script"
EOF
)

job_id=$(printf '%s\n' "$sbatch_script" | sbatch --parsable)

cat << EOF
submitted nextflow head process as job $job_id
   portal:    $mode
   partition: $queue   walltime: $time_limit   cpus: $cpus   mem: $mem
   log:       $log_folder/$job_id.out
   status:    scripts/status-pipeline.sh $job_id
EOF
