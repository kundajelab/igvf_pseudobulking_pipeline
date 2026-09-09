include { dotenv } from 'plugin/nf-dotenv'
include { oomMemoryOf ; oomMaxRetriesOf } from './retry.nf'

process SPLIT_FRAGMENTS {
    // -in environments with high contention, only use 5 cpus (any more will actually be *SLOWER*)
    // -in environments with low contention, use 8 cpus (finishes in a reasonable amount of time)
    cpus params.split_fragments_parallel ? 8 : 5
    memory { oomMemoryOf(task, baseMem) }
    maxRetries { oomMaxRetriesOf(task, baseMem) }
    conda "environments/PSEUDOBULK.yaml"
    container "${dotenv('PSEUDOBULK_IMAGE')}"

    input:
        tuple path(fragments_file), val(num_unique_barcodes), val(num_kept_fragments)
        path(metadata_file)
        path(chrom_sizes)
        path(tss_tsv)
    output:
        path("${local_output_folder}/separated_pseudorep1/*.tsv"), optional: true, emit: pseudorep_1
        path("${local_output_folder}/separated_pseudorep2/*.tsv"), optional: true, emit: pseudorep_2
        path("${local_output_folder}/separated_pseudorepT/*.tsv"), optional: true, emit: pseudorep_t
        path("${local_output_folder}/separated_fragments/*.tsv"), optional: true, emit: separated_fragments
        tuple val(analysis_set_accession), path("${local_output_folder}/atac_qc_reports/${analysis_set_accession}.tsv"), emit: atac_qc_reports
        tuple val(analysis_set_accession), path("${local_output_folder}/atac_qc_reports/${analysis_set_accession}_tss_matrix.npz"), emit: tss_matrix
        path("${local_output_folder}/atac_qc_reports/*"), emit: atac_qc_files

    script:
    analysis_set_accession = fragments_file.getSimpleName()
    local_output_folder = "output"
    // Two things scale with the input, both counted by COUNT_BARCODES:
    // - QC data is held in memory for every barcode (dominated by transcription start sites array)
    // - every pseudobulked fragment is held as a Fragment object, and they are all live at once
    //   because no pseudobulk is written out until all of them have been collected
    // NOTE: 200.B per fragment comes from measuring retained heap for the real class: ~197 bytes for
    // the slotted object, its uncached np.int32 start and end, and its share of the contig and
    // barcode strings. Fragment.from_line interns those two strings so that all the fragments in a
    // file share one object per distinct value; without that it measures ~346 bytes instead.
    baseMem = 1.GB + 12.KB * num_unique_barcodes + 200.B * num_kept_fragments
    """
    1>&2 echo "Requested ${task.memory} to process ${num_unique_barcodes} unique barcodes and ${num_kept_fragments} fragments."
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
