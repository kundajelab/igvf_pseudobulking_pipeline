include { dotenv } from 'plugin/nf-dotenv'
include { oomMemoryOf ; oomMaxRetriesOf } from './retry.nf'

process PSEUDOBULK_RNA {
    cpus { Math.min(4, h5ads.size()) }
    memory { oomMemoryOf(task, baseMem) }
    maxRetries { oomMaxRetriesOf(task, baseMem) }
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
    publishDir "${params.workspace}/${params.principal_analysis.replace(",", "-")}/output",
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
    // Memory model, fitted to the measured peak anonymous memory of the whole process tree over 11
    // runs spanning 1-8 workers and 166 MB - 7.1 GB of input:
    //
    //     peak = PARENT_OVERHEAD + RATIO * totalSize + workers * PER_WORKER
    //
    // Each worker holds one whole h5ad at a time, so the data-proportional term is the sum of the
    // sizes of the largest `workers` inputs rather than a multiple of the largest one: the input
    // sizes are badly skewed, and one big h5ad next to three small ones needs far less than four
    // big ones. Each worker is its own interpreter (the pool uses processes, not threads) holding
    // its own gene reference and slimmed metadata, measured at ~460 MB.
    totalSize = h5ads.collect { h5ad -> h5ad.size() }
        .sort{ s -> -s }.take(task.cpus).sum()
    baseMem = 1.MB * ((1024 + 7.5 * totalSize / 2 ** 20 + task.cpus * 512) as long)
    """
    export PYTHON_GIL=1
    pseudobulk pseudobulk-rna \
        --input input_dir \
        --output-dir . \
        --metadata-loc "${metadata_file}" \
        --gene-info "${gene_info}" \
        --num-workers ${task.cpus}
    """
}
