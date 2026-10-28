include { dotenv } from 'plugin/nf-dotenv'
include { oomMemoryOf ; oomMaxRetriesOf } from './retry.nf'

// Call peaks on the three pseudoreps of one pseudobulk in a single job, running the three side by
// side. One job rather than three keeps the many short MACS3 runs from running into the
// scheduler's limit on jobs submitted per hour, and running them in parallel keeps the job no
// longer than its slowest pseudorep (T, which holds the insertions of both 1 and 2).
// For readable logs, every line a pseudorep prints is prefixed with "[rep N]" as it is written, and
// the unprefixed output of each is also saved to macs3_rep_N.log in the work directory.
process MACS3 {
    cpus 3
    conda "environments/CALL_PEAKS.yaml"
    container "${dotenv('CALL_PEAKS_IMAGE')}"
    cache 'deep'
    memory { oomMemoryOf(task, baseMem) }
    maxRetries { oomMaxRetriesOf(task, baseMem) }

    input:
        tuple val(pseudobulk_id),
            path(pseudorep_1, arity: "1"),
            path(pseudorep_2, arity: "1"),
            path(pseudorep_t, arity: "1")
        path(chr_order)
        val(species)
    output:
        tuple val(pseudobulk_id),
            path("${pseudobulk_id}.1_peaks_top.narrowPeak"),
            path("${pseudobulk_id}.1_ppois_clipped.bdg"),
            path("${pseudobulk_id}.2_peaks_top.narrowPeak"),
            path("${pseudobulk_id}.2_ppois_clipped.bdg"),
            path("${pseudobulk_id}.t_peaks_top.narrowPeak"),
            path("${pseudobulk_id}.t_ppois_clipped.bdg"),
            emit: output

    script:
    // Heuristic on memory needed by one pseudorep: 2 GB + 6 * size of its fragments bed. The three
    // run at the same time, each in its own macs3 process, so they need the sum.
    baseMem = 1.MB * (
        [pseudorep_1, pseudorep_2, pseudorep_t]
            .collect { pseudorep -> 2048 + 6.0 * pseudorep.size() / 2 ** 20 }
            .sum() as long
    )
    // MACS3 effective genome size shortcut for each supported species
    genome_size = [human: "hs", mouse: "mm"][species]
    if (genome_size == null) {
        error "MACS3: no effective genome size for species '${species}'"
    }
    // Cap sort's buffer at a third of the task's memory, less some for the rest of the pipeline, as
    // the three pseudoreps' sorts can run at once. By default sort sizes it from the node's physical
    // memory, which can exceed the job's limit.
    sort_buffer_size = "${Math.max(256, task.memory.toMega().intdiv(3) - 512)}M"
    """
    # avoid running out of /tmp space in cluster/cloud environment
    # NOTE: one temporary folder for the whole job, removed by this shell's trap. A trap set inside
    # a pseudorep's background pipeline is not reliably run when it exits.
    temp_dir=\$(mktemp -d -p .)
    trap 'rm -rf "\$temp_dir"' EXIT

    # Call peaks on one pseudorep, then make its top peaks and clipped p-value bedgraph.
    function call_pseudorep_peaks {
        local -r rep="\$1"
        local -r fragments_tsv="\$2"
        local -r base="${pseudobulk_id}.\$rep"

        macs3 callpeak \
            -t "\$fragments_tsv" \
            -f BED \
            -n "\$base" \
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
        local -r rep_temp_dir="\$temp_dir/\$rep"
        mkdir "\$rep_temp_dir"

        # prevent sigpipe by sorting ascending by -log10(pvalue) and taking the last (strongest) values
        sort \
            --temporary-directory="\$rep_temp_dir" \
            --buffer-size="${sort_buffer_size}" \
            -k 8g,8g \
            "\${base}_peaks.narrowPeak" \
        | tail -n "${params.num_top_peaks}" \
        | sort-bed.sh --buffer-size "${sort_buffer_size}" "${chr_order}" \
        > "\${base}_peaks_top.narrowPeak"

        1>&2 echo "Making p-value bedgraphs"
        macs3 bdgcmp \
            -m ppois \
            -t "\${base}_treat_pileup.bdg" \
            -c "\${base}_control_lambda.bdg" \
            -o "\${base}_ppois.bdg"

        1>&2 echo "Clipping bedgraphs"
        bedClip \
            "\${base}_ppois.bdg" \
            "${chr_order}" \
            "\${base}_ppois_clipped.bdg"
    }

    # Start one pseudorep in the background. Its output (stdout, stderr and the shell trace) goes
    # to macs3_rep_<rep>.log, and to this job's stderr with each line prefixed by "[rep <rep>]".
    # NOTE: pipefail in the subshell makes its exit status that of call_pseudorep_peaks rather than
    # of sed, so a failure, such as being killed for running out of memory, is not lost.
    # NOTE: indexed arrays rather than associative ones, which bash 3 (e.g. macOS) does not have.
    reps=(1 2 t)
    pseudorep_tsvs=("${pseudorep_1}" "${pseudorep_2}" "${pseudorep_t}")
    pseudorep_pids=()
    for index in "\${!reps[@]}"; do
        rep="\${reps[\$index]}"
        (
            set -o pipefail
            call_pseudorep_peaks "\$rep" "\${pseudorep_tsvs[\$index]}" 2>&1 \
            | tee "macs3_rep_\$rep.log" \
            | sed -u "s/^/[rep \$rep] /" 1>&2
        ) &
        pseudorep_pids+=("\$!")
    done

    # Wait for all three, and fail with the exit status of the first that failed, so that an
    # out-of-memory kill (137) is still retried with more memory.
    failed_reps=()
    exit_status=0
    for index in "\${!reps[@]}"; do
        rep="\${reps[\$index]}"
        if wait "\${pseudorep_pids[\$index]}"; then
            1>&2 echo "[rep \$rep] finished"
        else
            rep_status=\$?
            failed_reps+=("\$rep")
            if [[ "\$exit_status" -eq 0 ]]; then
                exit_status="\$rep_status"
            fi
        fi
    done
    if [[ "\${#failed_reps[@]}" -gt 0 ]]; then
        1>&2 echo "MACS3 failed for pseudorep(s) \${failed_reps[*]} (first exit status \$exit_status):" \
            "see their [rep N] lines above, or macs3_rep_N.log in the work directory"
        exit "\$exit_status"
    fi
    """
}
