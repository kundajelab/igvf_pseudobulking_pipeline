include { dotenv } from 'plugin/nf-dotenv'

process SUMMARIZE_PSEUDOBULK_QC {
    cpus 1
    memory '8 GB'
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
        path("atac_qc_dir/*")
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
