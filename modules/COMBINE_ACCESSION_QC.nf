include { dotenv } from 'plugin/nf-dotenv'
include { oomMemoryOf ; oomMaxRetriesOf } from './retry.nf'

// Create a combined QC report for each intermediate analysis set accession
process COMBINE_ACCESSION_QC {
    cpus {
        def num_atac_qc_files = atac_qc_files.sum { atac_qc -> atac_qc.getExtension() == "npz" ? 0 : 1 } ?: 0
        def num_rna_qc_files = rna_qc_files.size()
        Math.min(
            16,
            Math.max(
                1,
                Math.max(num_atac_qc_files, num_rna_qc_files)
            )
        )
    }
    memory { oomMemoryOf(task, baseMem) }
    maxRetries { oomMaxRetriesOf(task, baseMem) }
    conda "environments/PSEUDOBULK.yaml"
    container "${dotenv('PSEUDOBULK_IMAGE')}"
    publishDir "${params.workspace}/${params.principal_analysis.replace(",", "-")}/output/analysis_accession_qc_reports",
        pattern: "*_per_cell_qc.tsv.gz",
        mode: params.publish_mode
    cache false

    input:
        path(atac_qc_files, name: "atac_qc_dir/*", arity: "0..*")
        path(rna_qc_files, name: "rna_qc_dir/*", arity: "0..*")
        path(metadata_file)
    output:
        path("*_per_cell_qc.tsv.gz"), emit: accession_qcs

    script:
    // NOTE: supply a zero initial value to each sum, because summing an empty list of files (which
    // happens e.g. when there were no fragments files, so no ATAC QC) otherwise yields null
    // We will hold a whole h5ad in path for each CPU, so take the
    // N largest h5ad files where N=number of cpus. Start sum at 0 in case RNA is missing.
    rna_qc_size = rna_qc_files.collect { rna_qc -> rna_qc.size() }
        .sort{s -> -s}.take(task.cpus).sum() ?: 0
    // similar for ATAC files
    atac_qc_size = atac_qc_files.collect { atac_qc -> atac_qc.getExtension() == "npz" ? 0 : atac_qc.size() }
        .sort{s -> -s}.take(task.cpus).sum() ?: 0

    // heuristic: we add pad with 3 GB baseline and estimate using around 3.5x space due to merging and duplication
    totalSize = rna_qc_size + atac_qc_size
    baseMem = 3.GB + 1.GB * Math.round(totalSize * 3.5 / 2 ** 30)
    """
    export PYTHON_GIL=1
    pseudobulk combine-accession-qc \
        --metadata-loc "${metadata_file}" \
        --atac-qc-dir "atac_qc_dir" \
        --rna-qc-dir "rna_qc_dir" \
        --output-dir . \
        --num-workers ${task.cpus}
    """
}
