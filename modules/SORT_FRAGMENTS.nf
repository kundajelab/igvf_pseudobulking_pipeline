include { dotenv } from 'plugin/nf-dotenv'

process SORT_FRAGMENTS {
    cpus 2
    memory '8 GB'
    conda "environments/CALL_PEAKS.yaml"
    container "${dotenv('CALL_PEAKS_IMAGE')}"
    publishDir "${params.workspace}/${params.principal_analysis.replace(",", "-")}/output/pseudobulks/${pseudobulk_id}",
        saveAs: { _file_name -> publish ? "fragments.tsv.gz" : null },
        mode: params.publish_mode

    input:
        tuple val(pseudobulk_id), path(fragments_tsvs)
        path(chrom_sizes)
        val(publish)
    output:
        tuple val(pseudobulk_id), path(sorted_fragments_tsv), emit: sorted_fragments_tsv

    script:
    sorted_fragments_tsv = "${fragments_tsvs[0].getBaseName(2)}.sorted.tsv.gz"
    """
    # sort the concatenated fragment TSVs using bin/sort-bed.sh, then bgzip
    sort-bed.sh "${chrom_sizes}" "${fragments_tsvs.join('" "')}" \
    | bgzip -@ ${task.cpus} -o "${sorted_fragments_tsv}"
    """
}
