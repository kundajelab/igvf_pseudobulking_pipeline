include { dotenv } from 'plugin/nf-dotenv'

// Create a combined QC report for each intermediate analysis set accession
process COMBINE_ACCESSION_QC {
    cpus { Math.min(16, Math.max(num_atac_qc_files, rna_qc_files.size())) }
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
    publishDir "${params.workspace}/${params.principal_analysis.replace(",", "-")}/output/analysis_accession_qc_reports",
        pattern: "*_per_cell_qc.tsv.gz",
        mode: params.publish_mode

    input:
        path(atac_qc_files, name: "atac_qc_dir/*", arity: "1..*")
        path(rna_qc_files, name: "rna_qc_dir/*", arity: "1..*")
        path(metadata_file)
    output:
        path("*_per_cell_qc.tsv.gz"), emit: accession_qcs

    script:
    rna_qc_size = rna_qc_files.sum { rna_qc -> rna_qc.size() }
    num_atac_qc_files = atac_qc_files.sum { atac_qc -> atac_qc.getExtension() == ".npz" ? 0 : 1 }
    atac_qc_size = atac_qc_files.sum { atac_qc -> atac_qc.getExtension() == ".npz" ? 0 : atac_qc.size() }
    totalSize = rna_qc_size + atac_qc_size
    baseMem = 1.GB + 1.GB * Math.round(totalSize * 1.25 / 2 ** 30)
    """
    export PYTHON_GIL=0
    pseudobulk combine-accession-qc \
        --metadata-loc "${metadata_file}" \
        --atac-qc-dir "atac_qc_dir" \
        --rna-qc-dir "rna_qc_dir" \
        --output-dir . \
        --num-workers ${task.cpus}
    """
}
