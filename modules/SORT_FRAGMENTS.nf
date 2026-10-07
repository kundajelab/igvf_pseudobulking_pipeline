include { dotenv } from 'plugin/nf-dotenv'
include { oomMemoryOf ; oomMaxRetriesOf } from './retry.nf'

process SORT_FRAGMENTS {
    cpus 2
    memory { oomMemoryOf(task, baseMem) }
    maxRetries { oomMaxRetriesOf(task, baseMem) }
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
    sorted_fragments_bigbed = "${fragments_tsvs[0].getBaseName(2)}.fragments.bb"
    // Memory model: sorting in memory peaks at 3.4x the uncompressed input (measured over 52 tasks of
    // all four kinds, 0.01 - 1.3 GB of input). Sorts that spilled to temp files once their input
    // outgrew the buffer were killed at 8 GB with a 4 GB buffer, so request enough to sort the
    // whole input in memory, and give nearly all of it to sort's buffer.
    inputSize = fragments_tsvs.collect { tsv -> tsv.size() }.sum()
    baseMem = 1.GB + 1.MB * ((3.5 * inputSize / 2 ** 20) as long)
    // Leave 1 GB for awk, cut, bgzip and the shell. By default sort sizes its buffer from the node's
    // physical memory, which can exceed the job's limit.
    sort_buffer_size = "${Math.max(256, task.memory.toMega() - 1024)}M"
    """
    # sort the concatenated fragment TSVs using bin/sort-bed.sh, then bgzip
    sort-bed.sh --buffer-size "${sort_buffer_size}" "${chrom_sizes}" "${fragments_tsvs.join('" "')}" \
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
