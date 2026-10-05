include { dotenv } from 'plugin/nf-dotenv'
include { oomMemoryOf ; oomMaxRetriesOf } from './retry.nf'

process MACS3 {
    cpus 1
    conda "environments/CALL_PEAKS.yaml"
    container "${dotenv('CALL_PEAKS_IMAGE')}"
    cache 'deep'
    memory { oomMemoryOf(task, baseMem) }
    maxRetries { oomMaxRetriesOf(task, baseMem) }

    input:
        tuple val(pseudobulk_id), path(fragments_tsv, arity: "1")
        path(chr_order)
        val(species)
        val(suffix)
    output:
        tuple val(pseudobulk_id), path(top_peak_calls), path(clipped_ppois), emit: output

    script:
    fragmentsSize = fragments_tsv.size()
    // heuristic on memory needed: 2 + 6 * size of fragments bed, scaled by task attempt
    baseMem = 2.GB + 1.GB * (6.0 * fragmentsSize / 2 ** 30)
    // MACS3 effective genome size shortcut for each supported species
    genome_size = [human: "hs", mouse: "mm"][species]
    if (genome_size == null) {
        error "MACS3: no effective genome size for species '${species}'"
    }
    // Cap sort's buffer at half the task's memory, leaving room for the rest of the pipeline. By
    // default sort sizes it from the node's physical memory, which can exceed the job's limit.
    sort_buffer_size = "${task.memory.toMega().intdiv(2)}M"
    base = "${pseudobulk_id}.${suffix}"
    raw_peak_calls = "${base}_peaks.narrowPeak"
    top_peak_calls = "${base}_peaks_top.narrowPeak"
    treatment = "${base}_treat_pileup.bdg"
    control = "${base}_control_lambda.bdg"
    ppois = "${base}_ppois.bdg"
    clipped_ppois = "${base}_ppois_clipped.bdg"
    """
    macs3 callpeak \
        -t "${fragments_tsv}" \
        -f BED \
        -n "${base}" \
        -g "${genome_size}" \
        --outdir . \
        -p 0.01 \
        --shift -75 \
        --extsize 150 \
        --nomodel \
        -B \
        --SPMR \
        --keep-dup all \
        --call-summits

    1>&2 echo "Subsetting to top peaks"
    # avoid running out of /tmp space in cluster/cloud environment
    temp_dir=\$(mktemp -d -p .)
    trap 'rm -rf "\$temp_dir"' EXIT

    # prevent sigpipe by sorting ascending by -log10(pvalue) and taking the last (strongest) values
    sort \
        --temporary-directory="\$temp_dir" \
        --buffer-size="${sort_buffer_size}" \
        -k 8g,8g \
        "${raw_peak_calls}" \
    | tail -n "${params.num_top_peaks}" \
    | sort-bed.sh --buffer-size "${sort_buffer_size}" "${chr_order}" \
    > "${top_peak_calls}"

    1>&2 echo "Making p-value bedgraphs"
    macs3 bdgcmp \
        -m ppois \
        -t "${treatment}" \
        -c "${control}" \
        -o "${ppois}"

    1>&2 echo "clipping bedgraphs"
    bedClip \
        "${ppois}" \
        "${chr_order}" \
        "${clipped_ppois}"
    """


}
