include { dotenv } from 'plugin/nf-dotenv'
include { oomMemoryOf ; oomMaxRetriesOf } from './retry.nf'

// Sort all four fragment sets of one pseudobulk in a single job: the separated fragments that are
// published, and the three pseudoreps that peaks are called on. Each set gathers that pseudobulk's
// files from every analysis set. One job rather than four keeps the many short sorts from running
// into the scheduler's limit on jobs submitted per hour.
process SORT_FRAGMENTS {
    cpus 2
    memory { oomMemoryOf(task, baseMem) }
    maxRetries { oomMaxRetriesOf(task, baseMem) }
    conda "environments/CALL_PEAKS.yaml"
    container "${dotenv('CALL_PEAKS_IMAGE')}"
    publishDir "${params.workspace}/${params.principal_analysis.replace(",", "-")}/output/pseudobulks/${pseudobulk_id}",
        pattern: "${sorted_fragments_tsv}",
        saveAs: { _file_name -> "fragments.tsv.gz" },
        mode: params.publish_mode
    publishDir "${params.workspace}/${params.principal_analysis.replace(",", "-")}/output/pseudobulks/${pseudobulk_id}",
        pattern: "${sorted_fragments_bigbed}",
        saveAs: { _file_name -> "fragments.bb" },
        mode: params.publish_mode

    input:
        // NOTE: arity "1..*" keeps each a list even for a pseudobulk from a single analysis set.
        // Otherwise that one file is bound as a bare path, which groovy iterates as its name parts.
        // The four sets stage without clashing, as split-fragments names them
        // <pseudobulk>.<analysis set>.tsv, and .1.tsv, .2.tsv and .t.tsv for the pseudoreps.
        tuple val(pseudobulk_id),
            path(fragments_tsvs, arity: "1..*"),
            path(pseudorep_1_tsvs, arity: "1..*"),
            path(pseudorep_2_tsvs, arity: "1..*"),
            path(pseudorep_t_tsvs, arity: "1..*")
        path(chrom_sizes)
        path(fragments_auto_sql)
    output:
        tuple val(pseudobulk_id), path(sorted_fragments_tsv), emit: sorted_fragments_tsv
        path(sorted_fragments_bigbed), emit: fragments_bigbed
        tuple val(pseudobulk_id), path(sorted_pseudorep_1), emit: pseudorep_1
        tuple val(pseudobulk_id), path(sorted_pseudorep_2), emit: pseudorep_2
        tuple val(pseudobulk_id), path(sorted_pseudorep_t), emit: pseudorep_t

    script:
    sorted_fragments_tsv = "${pseudobulk_id}.sorted.tsv.gz"
    sorted_fragments_bigbed = "${pseudobulk_id}.fragments.bb"
    sorted_pseudorep_1 = "${pseudobulk_id}.1.sorted.tsv.gz"
    sorted_pseudorep_2 = "${pseudobulk_id}.2.sorted.tsv.gz"
    sorted_pseudorep_t = "${pseudobulk_id}.t.sorted.tsv.gz"
    // Memory model: sorting in memory peaks at 3.4x the uncompressed input (measured over 52 tasks of
    // all four kinds, 0.01 - 1.3 GB of input). Sorts that spilled to temp files once their input
    // outgrew the buffer were killed at 8 GB with a 4 GB buffer, so request enough to sort the
    // whole input in memory, and give nearly all of it to sort's buffer. The four sorts run one
    // after another, so only the largest set has to fit.
    inputSize = [fragments_tsvs, pseudorep_1_tsvs, pseudorep_2_tsvs, pseudorep_t_tsvs]
        .collect { tsvs -> tsvs.collect { tsv -> tsv.size() }.sum() }
        .max()
    baseMem = 1.GB + 1.MB * ((3.5 * inputSize / 2 ** 20) as long)
    // Leave 1 GB for awk, cut, bgzip and the shell. By default sort sizes its buffer from the node's
    // physical memory, which can exceed the job's limit.
    sort_buffer_size = "${Math.max(256, task.memory.toMega() - 1024)}M"
    """
    # sort each set of concatenated fragment TSVs using bin/sort-bed.sh, then bgzip
    function sort_fragments {
        local -r sorted="\$1"
        shift 1
        sort-bed.sh --buffer-size "${sort_buffer_size}" "${chrom_sizes}" "\$@" \
        | bgzip -@ ${task.cpus} -o "\$sorted"
    }

    sort_fragments "${sorted_fragments_tsv}" "${fragments_tsvs.join('" "')}"
    # make bigbed version of the fragments BED, which is published along with it
    bedToBigBed \
        -type=bed3+2 \
        -as="${fragments_auto_sql}" \
        "${sorted_fragments_tsv}" \
        "${chrom_sizes}" \
        "${sorted_fragments_bigbed}"

    # NOTE: the pseudoreps are only used for peak calling, and macs3 can handle unsorted input.
    # However sorting is fast and makes the order consistent, so sort them anyway.
    sort_fragments "${sorted_pseudorep_1}" "${pseudorep_1_tsvs.join('" "')}"
    sort_fragments "${sorted_pseudorep_2}" "${pseudorep_2_tsvs.join('" "')}"
    sort_fragments "${sorted_pseudorep_t}" "${pseudorep_t_tsvs.join('" "')}"
    """
}
