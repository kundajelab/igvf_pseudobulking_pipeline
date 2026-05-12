include { dotenv } from 'plugin/nf-dotenv'

process SPLIT_FRAGMENTS {
    // -in environments with high contention, only use 5 cpus (any more will actually be *SLOWER*)
    // -in environments with low contention, use 8 cpus (finishes in a reasonable amount of time)
    cpus params.split_fragments_parallel ? 8 : 5
    memory "${4.GB * task.cpus}"
    conda "environments/PSEUDOBULK.yaml"
    container "${dotenv('PSEUDOBULK_IMAGE')}"
    cache "deep"

    input:
        path(fragments_file)
        path(metadata_file)
        path(chrom_sizes)
        path(tss_tsv)
    output:
        path("${local_output_folder}/separated_pseudorep1/*.tsv"), emit: pseudorep_1
        path("${local_output_folder}/separated_pseudorep2/*.tsv"), emit: pseudorep_2
        path("${local_output_folder}/separated_pseudorepT/*.tsv"), emit: pseudorep_t
        path("${local_output_folder}/separated_fragments/*.tsv"), emit: separated_fragments
        tuple val(analysis_set_accession), path("${local_output_folder}/atac_qc_reports/${analysis_set_accession}.tsv"), emit: atac_qc_reports
        tuple val(analysis_set_accession), path("${local_output_folder}/atac_qc_reports/${analysis_set_accession}_tss_matrix.npz"), emit: tss_matrix
        path("${local_output_folder}/atac_qc_reports/*"), emit: atac_qc_files

    script:
    analysis_set_accession = fragments_file.getSimpleName()
    local_output_folder = "output"
    """
    export PYTHON_GIL=0
    export PYTHON_JIT=1
    pseudobulk split-fragments \
        --fragments-file "${fragments_file}" \
        --output-dir "${local_output_folder}" \
        --metadata-loc "${metadata_file}" \
        --chrom-sizes "${chrom_sizes}" \
        --tss-tsv "${tss_tsv}" \
        --num-threads ${task.cpus}
    """
}
