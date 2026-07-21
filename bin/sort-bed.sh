#!/usr/bin/env bash
set -euo pipefail

# Sort BED file by chromosome (as specified by genome file), pos, and end
# 1. read in the chromosome order, assigning chrom-index based on order in the file
# 2. read from the BED file, and insert a new 1st field that is index into the chrom order
# 3. sort numerically based on chrom-index, start, end, and finally sort whatever is left
#    alphabetically
# 4. remove the chrom index

# Usage: sort-bed.sh <genome_order> [<bed>]
# If bed is not specified (or "-"), reads from stdin
genome_order=$1
shift 1

# avoid running out of /tmp space in cluster/cloud environment by making temp in working folder
temp_dir=$(mktemp -d -p .)
trap 'rm -rf "$temp_dir"' EXIT

# function to transparently read from stdin or compressed/uncompressed input BED
function read_bed {
    local -r bed="$1"
    if [[ "$bed" == "-" ]]; then
        cat
    elif [[ "$bed" =~ .gz ]]; then
        zcat < "$bed"
    else
        cat < "$bed"
    fi
}

function read_beds {
    if [[ "$#" == 0 ]]; then
        read_bed "-"
    else
        for bed in "${@}"; do
            read_bed "$bed"
        done
    fi
}

# 1. read from genome_order to build an index of chromosome order
# 2. then read from BED and add an initial column that is chromosome order
# 3. sort by chromosome index, start, end, then any remaining columns
# 4. cut away the chromosome index to yield sorted bed rows
awk -v OFS='\t' '
    FNR == 1 { ++file_num }
    file_num == 1 { idx[$1] = FNR }
    file_num == 2 { print idx[$1], $0 }
' "$genome_order" <(read_beds "${@}") \
| sort --temporary-directory="$temp_dir" -k1,1n -k3,3n -k4,4n -k5 \
| cut -f2-
