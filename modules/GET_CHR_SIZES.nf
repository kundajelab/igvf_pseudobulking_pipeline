include { dotenv } from 'plugin/nf-dotenv'

process GET_CHR_SIZES {
    cpus 2
    memory '8 GB'
    conda "environments/CALL_PEAKS.yaml"
    container "${dotenv('CALL_PEAKS_IMAGE')}"

    input:
        path(fasta)
    output:
        tuple val(fasta.simpleName), path(chr_sizes), emit: chr_sizes

    script:
    fasta_name = fasta.simpleName
    chr_sizes = "${fasta_name}.chr_sizes.tsv"
    """
    bgzip -cd --threads ${task.cpus} "${fasta}" \
    | samtools faidx --threads ${task.cpus} --fai-idx /dev/stdout - \
    | cut -f 1,2 \
    > "${chr_sizes}"
    """
}
