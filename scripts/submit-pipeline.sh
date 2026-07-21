#!/usr/bin/env bash
set -euo pipefail
# WARNING: this doesn't work because sbatch scripts can't run on dev queue. Need to rethink and
# probably use srun instead of sbatch

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
repo_dir=$(cd -- "$script_dir/.." &> /dev/null && pwd)

metadata_file=$(readlink -f "$1")
queue=${2:-dev}

temp_dir=$(mktemp -d)
trap 'rm -rf "$temp_dir"' EXIT
sbatch_script="$temp_dir/igvf-pseudobulk.sbatch"
name="igvf-pseudobulking/$(basename "$metadata_file")"
output="$HOME/logs/${name%.gz}"
output="${output%.tsv}.out"
mkdir -p "$(dirname "$output")"

cat << EOF > "$sbatch_script"
#!/usr/bin/env bash
#SBATCH --job-name=$name
#SBATCH --output=$output
#SBATCH --time=12:00:00
#SBATCH --partition=$queue
#SBATCH --cpus-per-task=1 
#SBATCH --mem=20GB
#SBATCH --chdir="$repo_dir"
hostname
pixi run pipeline "$metadata_file"
EOF

1>&2 echo "$sbatch_script:"
1>&2 cat "$sbatch_script"
sbatch "$sbatch_script"
