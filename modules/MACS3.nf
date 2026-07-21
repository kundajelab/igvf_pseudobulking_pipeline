include { dotenv } from 'plugin/nf-dotenv'

process MACS3 {
    cpus 1
    conda "environments/CALL_PEAKS.yaml"
    container "${dotenv('CALL_PEAKS_IMAGE')}"
    cache 'deep'
    memory {
        if (task.previousTrace) {
            def wasOom = task.exitStatus in 137..140
            wasOom ? task.previousTrace.memory + baseMem : baseMem
        } else {
            baseMem
        }
    }
    maxRetries {
        def wasOom = task.exitStatus in 137..140
        def previousOomCount = (task.previousTrace.memory / baseMem).round() as Integer
        task.executor == 'local' || (
            previousOomCount >= params.max_oom_retries && wasOom
        ) ?
        params.oom_max_retries :
        params.preemptible_max_retries + params.oom_max_retries
    }

    input:
        tuple val(pseudobulk_id), path(fragments_tsv, arity: "1")
        path(chr_order)
        val(suffix)
    output:
        tuple val(pseudobulk_id), path(top_peak_calls), path(clipped_ppois), emit: output

    script:
    fragmentsSize = fragments_tsv.size()
    // heuristic on memory needed: 2 + 6 * size of fragments bed, scaled by task attempt
    baseMem = 2.GB + 1.GB * (6.0 * fragmentsSize / 2 ** 30)
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
        -g hs \
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

    # prevent sigpipe by reversing the sort and taking the last values
    sort \
        --reverse \
        --temporary-directory="\$temp_dir" \
        -k 8gr,8gr  \
        "${raw_peak_calls}" \
    | tail -n "${params.num_top_peaks}" \
    | sort-bed.sh "${chr_order}" \
    > "${top_peak_calls}"

    1>&2 echo "Making p-value bedgraphs"
    macs3 bdgcmp \
        -m ppois \
        -t "${treatment}" \
        -c "${control}" \
        -o "${ppois}"
    
    1>2 echo "clipping bedgraphs"
    bedClip \
        "${ppois}" \
        "${chr_order}" \
        "${clipped_ppois}"
    """


}
