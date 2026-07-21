include { dotenv } from 'plugin/nf-dotenv'

process PSEUDOBULK_RNA {
    cpus 4
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
    conda "environments/PSEUDOBULK.yaml"
    container "${dotenv('PSEUDOBULK_IMAGE')}"
    cache "deep"
    publishDir "${params.workspace}/${params.principal_analysis.replace(",", "-")}/output/pseudobulks",
        pattern: "pseudobulks/*",
        saveAs: { file_name ->
            def pseudobulk_id = file_name.tokenize("/")[-1].tokenize(".")[0]
            file_name.endsWith(".h5ad") ?
                "${pseudobulk_id}/rna_counts_mtx.h5ad" :
                "${pseudobulk_id}/pseudobulk_expression.tsv.gz"
        },
        mode: params.publish_mode
    publishDir "${params.workspace}/output",
        pattern: "cell_name_to_annotation_mapping.tsv",
        mode: params.publish_mode

    input:
        path(h5ads, name: "input_dir/*", arity: "1..*")
        path(metadata_file)
        path(gene_info)
    output:
        path("rna_qc_reports/*.scRNA_all_cells_QC_metrics.tsv"), emit: all_cell_rna_qc_reports
        path("rna_qc_reports/*.pseudobulked_cell_QC_metrics.tsv"), emit: pseudobulked_rna_qc_reports
        path("pseudobulks/*.tsv.gz"), emit: pseudobulk_counts
        path("pseudobulks/*.h5ad"), emit: pseudobulk_h5ads
        path("cell_name_to_annotation_mapping.tsv"), emit: cell_name_to_annotation_mapping

    script:
    totalSize = h5ads.sum { h5ad -> h5ad.size() }
    baseMem = 1.GB + 1.GB * Math.round(totalSize * 1.25 / 2 ** 30)
    """
    export PYTHON_GIL=0
    pseudobulk pseudobulk-rna \
        --input input_dir \
        --output-dir . \
        --metadata-loc "${metadata_file}" \
        --gene-info "${gene_info}" \
        --num-workers ${task.cpus}
    """
}
