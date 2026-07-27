include { dotenv } from 'plugin/nf-dotenv'

process SUMMARIZE_PSEUDOBULK_QC {
    cpus 1
    // Peak RSS measured over 609 tasks was 1.1 GB typical, 1.5 GB worst case, so 8 GB was ~5x
    // oversubscribed and made these jobs needlessly hard to backfill on the owners queue. Grow on
    // retry so an unusually large pseudobulk still gets through instead of failing outright.
    memory { 2.GB * task.attempt }
    // This task should take around one minute. On rare occasions jobs can hang due to the workdir not having
    // files updated, and since there are many of these jobs they are more likely to be affected.
    // Time cap scales with attempt so a genuinely slow task gets room on retry.
    time { 15.min * task.attempt }
    conda "environments/PSEUDOBULK.yaml"
    container "${dotenv('PSEUDOBULK_IMAGE')}"
    publishDir "${params.workspace}/${params.principal_analysis.replace(",", "-")}/output/pseudobulks/${pseudobulk_id}",
        pattern: "${pseudobulk_qc_out}",
        saveAs: { _f -> "per_cell_qc.tsv.gz" },
        mode: params.publish_mode

    input:
        tuple val(pseudobulk_id),
            path(rna_qc),
            path(pseudobulk_counts),
            path(frip_per_cell),
            path(fragments_per_cell),
            path(fragments_in_peaks_per_cell)
        path("atac_qc_dir/*", arity: "0..*")
        path(metadata_file)
    output:
        path(pseudobulk_qc_out), emit: pseudobulk_qc_out
        path(qc_summary_out), emit: qc_summary_out

    script:
    pseudobulk_qc_out = "${pseudobulk_id}.per_cell_qc.tsv.gz"
    qc_summary_out = "${pseudobulk_id}.pseudobulk_qc.tsv"
    """
    export PYTHON_GIL=1
    pseudobulk summarize-pseudobulk-qc \
        --pseudobulk "${pseudobulk_id}" \
        --metadata-loc "${metadata_file}" \
        --atac-qc-dir "atac_qc_dir" \
        ${rna_qc.size() == 0 ? "" : "--rna-qc \"${rna_qc}\""} \
        ${pseudobulk_counts.size() == 0 ? "" : "--pseudobulk-counts \"${pseudobulk_counts}\""} \
        ${frip_per_cell.size() == 0 ? "" : "--frip-per-cell \"${frip_per_cell}\""} \
        ${fragments_per_cell.size() == 0 ? "" : "--fragments-per-cell \"${fragments_per_cell}\""} \
        ${fragments_in_peaks_per_cell.size() == 0 ? "" : "--fragments-in-peaks-per-cell \"${fragments_in_peaks_per_cell}\""} \
        --pseudobulk-qc-out "${pseudobulk_qc_out}" \
        --qc-summary-out "${qc_summary_out}"
    """
}
