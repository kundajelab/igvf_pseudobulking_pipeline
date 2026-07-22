include { dotenv } from 'plugin/nf-dotenv'

process SORT_FRAGMENTS {
    cpus 2
    memory '8 GB'
    conda "environments/CALL_PEAKS.yaml"
    container "${dotenv('CALL_PEAKS_IMAGE')}"
    publishDir "${params.workspace}/${params.principal_analysis.replace(",", "-")}/output/pseudobulks/${pseudobulk_id}",
        pattern: "${sorted_fragments_tsv}",
        saveAs: { _file_name -> publish ? "fragments.tsv.gz" : null },
        mode: params.publish_mode
    publishDir "${params.workspace}/${params.principal_analysis.replace(",", "-")}/output/pseudobulks/${pseudobulk_id}",
        pattern: "${sorted_fragments_bigbed}",
        saveAs: { _file_name -> publish ? "fragments.bb" : null },
        mode: params.publish_mode

    input:
        tuple val(pseudobulk_id), path(fragments_tsvs)
        path(chrom_sizes)
        path(fragments_auto_sql)
        val(publish)
    output:
        tuple val(pseudobulk_id), path(sorted_fragments_tsv), emit: sorted_fragments_tsv
        path(sorted_fragments_bigbed), optional: true, emit: fragments_bigbed

    script:
    sorted_fragments_tsv = "${fragments_tsvs[0].getBaseName(2)}.sorted.tsv.gz"
    sorted_fragments_bigbed = "${fragments_tsvs[0].getBaseName(2)}.bb"
    """
    # sort the concatenated fragment TSVs using bin/sort-bed.sh, then bgzip
    sort-bed.sh "${chrom_sizes}" "${fragments_tsvs.join('" "')}" \
    | bgzip -@ ${task.cpus} -o "${sorted_fragments_tsv}"

    if [[ "${publish}" == "true" ]]; then
        # make bigbed version of fragments BED
        bedToBigBed \
            -type=bed3+2 \
            -as="${fragments_auto_sql}" \
            "${sorted_fragments_tsv}" \
            "${chrom_sizes}" \
            "${sorted_fragments_bigbed}"
    fi
    """
}
