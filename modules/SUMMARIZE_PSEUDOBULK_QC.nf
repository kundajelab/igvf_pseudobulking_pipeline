include { dotenv } from 'plugin/nf-dotenv'

// Group a channel of per-pseudobulk tuples of the form
//   (pseudobulk_id, rna_qc, pseudobulk_counts, frip_per_cell, fragments_per_cell, fragments_in_peaks_per_cell)
// with an empty list standing in for each missing file, into batches of batch_size pseudobulks for
// SUMMARIZE_PSEUDOBULK_QC. Yields tuples of (manifest, files), where the manifest has one row of file
// names per pseudobulk ("" for missing files) and files are all the files of the batch.
// NOTE: the input channel emits in whatever order the upstream tasks finish, so sort by
// pseudobulk_id first to make the batches the same on every run; otherwise -resume would re-run
// most of them. This waits for every upstream task before the first batch is emitted.
def batchPseudobulkQcInputs(qc_inputs_ch, batch_size) {
    return qc_inputs_ch
        .toSortedList { a, b -> a[0] <=> b[0] }
        .flatMap()
        .buffer(size: batch_size, remainder: true)
        .map { batch ->
            def manifest = batch.collect { tup ->
                [tup[0]] + tup[1..5].collect { qc_file -> (qc_file instanceof List) ? "" : qc_file.name }
            }
            def files = batch.collectMany { tup -> tup[1..5].findAll { qc_file -> !(qc_file instanceof List) } }
            [manifest, files]
        }
}

// Summarize the QC for a batch of pseudobulks, one after another. The manifest has one row per
// pseudobulk of the form:
//   [pseudobulk_id, rna_qc, pseudobulk_counts, frip_per_cell, fragments_per_cell, fragments_in_peaks_per_cell]
// where each file is the name of a file in qc_inputs, or "" if that file is missing. Every input
// file name starts with its pseudobulk ID, so the files of a whole batch can be staged together.
// NOTE: missing files are written to the manifest TSV as "-" because tab is an IFS whitespace
// character, so "read" would collapse consecutive tabs and shift the later fields.
process SUMMARIZE_PSEUDOBULK_QC {
    cpus 1
    // Peak RSS measured over 609 tasks was 1.1 GB typical, 1.5 GB worst case, so 8 GB was ~5x
    // oversubscribed and made these jobs needlessly hard to backfill on the owners queue. Grow on
    // retry so an unusually large pseudobulk still gets through instead of failing outright.
    // The pseudobulks in a batch are summarized sequentially, so memory does not scale with batch size.
    memory { 2.GB * task.attempt }
    // Each pseudobulk should take around one minute. On rare occasions jobs can hang due to the workdir
    // not having files updated, and since there are many of these jobs they are more likely to be
    // affected; a hang now costs a retry of the whole batch. Time cap scales with attempt so a
    // genuinely slow task gets room on retry.
    time { (10.min + 2.min * manifest.size()) * task.attempt }
    conda "environments/PSEUDOBULK.yaml"
    container "${dotenv('PSEUDOBULK_IMAGE')}"
    publishDir "${params.workspace}/${params.principal_analysis.replace(",", "-")}/output/pseudobulks",
        pattern: "*.per_cell_qc.tsv.gz",
        saveAs: { file_name -> "${file_name.tokenize("/")[-1].tokenize(".")[0]}/per_cell_qc.tsv.gz" },
        mode: params.publish_mode

    input:
        tuple val(manifest), path(qc_inputs, name: "qc_inputs/*", arity: "0..*")
        path("atac_qc_dir/*", arity: "0..*")
        path(metadata_file)
    output:
        path "*.per_cell_qc.tsv.gz", arity: "1..*", emit: pseudobulk_qc_out
        path "*.pseudobulk_qc.tsv", arity: "1..*", emit: qc_summary_out

    script:
    manifest_tsv = manifest.collect { row -> row.collect { field -> field ?: "-" }.join("\t") }.join("\n")
    """
    export PYTHON_GIL=1
    cat > manifest.tsv << 'EOF'
${manifest_tsv}
EOF

    while IFS=\$'\\t' read -r pseudobulk_id rna_qc pseudobulk_counts frip_per_cell fragments_per_cell fragments_in_peaks_per_cell; do
        # NOTE: nextflow runs with "bash -u", under which bash < 4.4 treats expanding an empty array as
        # an unbound variable, hence the "+" guard where optional_args is expanded below.
        optional_args=()
        if [[ "\$rna_qc" != "-" ]]; then
            optional_args+=(--rna-qc "qc_inputs/\$rna_qc")
        fi
        if [[ "\$pseudobulk_counts" != "-" ]]; then
            optional_args+=(--pseudobulk-counts "qc_inputs/\$pseudobulk_counts")
        fi
        if [[ "\$frip_per_cell" != "-" ]]; then
            optional_args+=(--frip-per-cell "qc_inputs/\$frip_per_cell")
        fi
        if [[ "\$fragments_per_cell" != "-" ]]; then
            optional_args+=(--fragments-per-cell "qc_inputs/\$fragments_per_cell")
        fi
        if [[ "\$fragments_in_peaks_per_cell" != "-" ]]; then
            optional_args+=(--fragments-in-peaks-per-cell "qc_inputs/\$fragments_in_peaks_per_cell")
        fi
        1>&2 echo "Summarizing QC for \$pseudobulk_id"
        pseudobulk summarize-pseudobulk-qc \\
            --pseudobulk "\$pseudobulk_id" \\
            --metadata-loc "${metadata_file}" \\
            --atac-qc-dir "atac_qc_dir" \\
            \${optional_args[@]+"\${optional_args[@]}"} \\
            --pseudobulk-qc-out "\$pseudobulk_id.per_cell_qc.tsv.gz" \\
            --qc-summary-out "\$pseudobulk_id.pseudobulk_qc.tsv" \\
        < /dev/null
    done < manifest.tsv
    """

    stub:
    manifest_tsv = manifest.collect { row -> row.collect { field -> field ?: "-" }.join("\t") }.join("\n")
    """
    cat > manifest.tsv << 'EOF'
${manifest_tsv}
EOF

    while IFS=\$'\\t' read -r pseudobulk_id _rest; do
        touch "\$pseudobulk_id.per_cell_qc.tsv.gz" "\$pseudobulk_id.pseudobulk_qc.tsv"
    done < manifest.tsv
    """
}
