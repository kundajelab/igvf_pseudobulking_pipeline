include { dotenv } from 'plugin/nf-dotenv'

// Use aria2c to download all accession files in parallel
process DOWNLOAD_ACCESSION_FILES {
    secret 'IGVF_API_KEY'
    secret 'IGVF_SECRET_KEY'
    cpus 1
    memory '8 GB'
    conda "environments/DOWNLOAD_ACCESSION_FILES.yaml"
    container "${dotenv('DOWNLOAD_ACCESSION_FILES_IMAGE')}"
    cache "lenient"
    maxForks 1
    publishDir "${params.workspace}/${params.principal_analysis.replace(",", "-")}/raw_rna", pattern: "*.h5ad", mode: params.publish_mode
    publishDir "${params.workspace}/${params.principal_analysis.replace(",", "-")}/raw_fragments", pattern: "*.bed.gz", mode: params.publish_mode

    input:
        val(accessions)
    output:
        path("*.bed.gz"), emit: fragments_files
        path("*.h5ad"), emit: counts_matrix_files

    script:
    input_file = "aria2c_input.txt"
    """
    # get the download URLs (already redirected into S3) for files in each accession, and write to an aria2c input file
    for accession in ${accessions.join(' ')}; do
        get-url "\$accession" --output "${input_file}"
    done

    # download the files:
    # - get URLs and downloaded file names from the aria2c input_file
    # - split each download into 16 parts that can be downloaded independently
    # - download via at most 16 simultaneous connections per server
    # - download at most 8 files at once
    # - respect the system's open file limit
    aria2c \
        --input-file="${input_file}" \
        --max-connection-per-server=16 \
        --split=16 \
        --max-concurrent-downloads=8 \
        --rlimit-nofile \$(ulimit -n)
    """
}
