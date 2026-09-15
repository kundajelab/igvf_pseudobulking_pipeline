include { FIND_BAD_BEDS } from './modules/FIND_BAD_BEDS.nf'
include { DOWNLOAD_ACCESSION_FILES as DOWNLOAD_BEDS } from './modules/DOWNLOAD_ACCESSION_FILES.nf'
include { DOWNLOAD_ACCESSION_FILES as DOWNLOAD_FASTAS } from './modules/DOWNLOAD_ACCESSION_FILES.nf'
include { GET_CHR_SIZES } from './modules/GET_CHR_SIZES.nf'
include { FILTER_TO_BIG_BED } from './modules/FILTER_TO_BIG_BED.nf'
include { UPLOAD_FIXED_BEDS } from './modules/UPLOAD_FIXED_BEDS.nf'

workflow {
    // Find all the BED files that have problems
    FIND_BAD_BEDS(params.igvf_mode, params.max_pseudobulks, params.lab)

    // Download all the BED files that need to be fixed, download_batch_size accessions at a time
    // simultaneously. This limit:
    // 1) avoids exceeding the time limit on the redirected S3 paths
    // 2) allows downstream processing of the early batches while later batches are still
    //    downloading
    bed_accessions_ch = FIND_BAD_BEDS.out.bed_records
        .flatten()
        .map { bed -> "${bed.simpleName};${bed.simpleName}" }
        .toSortedList()
        .flatten()
    DOWNLOAD_BEDS(
        bed_accessions_ch.buffer(size: params.download_batch_size, remainder: true),
        params.igvf_mode
    )

    // Download all the FASTAs that are needed
    bed_references_ch = FIND_BAD_BEDS.out.ref_accession_files
        .flatten()
        .map { ref_id_file -> [ref_id_file.text.trim(), ref_id_file.simpleName] }
    ref_accessions_ch = bed_references_ch
        .map { ref_accession, _bed_accession -> ref_accession }
        .unique()
        .toSortedList()
        .flatten()
    DOWNLOAD_FASTAS(
        ref_accessions_ch
            .map { ref_a -> "${ref_a};${ref_a}" }
            .buffer(size: 16, remainder: true),
        params.igvf_mode
    )
    // Get the chr_sizes TSV for each reference, this is what the fixing Process needs.
    GET_CHR_SIZES(DOWNLOAD_FASTAS.out.fasta_files.flatten())


    // Combine all the downloaded BED files, chr_sizes and records JSONs into a channel that can be
    // used to fix the BEDs
    bed_records_ch = FIND_BAD_BEDS.out.bed_records
        .flatten()
        .map { bed -> [bed.simpleName, bed] }
    big_bed_records_ch = FIND_BAD_BEDS.out.big_bed_records
        .flatten()
        .map { bed -> [bed.simpleName, bed] }
    records_ch = bed_records_ch.join(big_bed_records_ch, remainder: true)
        .map { accession, bed_json, big_bed_json -> [accession, bed_json, big_bed_json] }
    fragments_files_ch = DOWNLOAD_BEDS.out.fragments_files
        .flatten()
        .map { bed -> [bed.simpleName, bed] }
    peaks_files_ch = DOWNLOAD_BEDS.out.peaks_files
        .flatten()
        .map { bed -> [bed.simpleName, bed] }
    // Several BEDs may share a reference. combine keeps every BED-to-FASTA match.
    bed_fastas_ch = bed_references_ch
        .combine(GET_CHR_SIZES.out.chr_sizes, by: 0)
        .map { _ref_accession, bed_accession, fasta -> [bed_accession, fasta] }

    // [bed_accession, bed_file, bed_json, big_bed_json (or []), chr_sizes]
    bed_ch = fragments_files_ch.concat(peaks_files_ch)
        .join(records_ch)
        .join(bed_fastas_ch)
    // We want to batch this rather than spinning up a task per BED
    bed_batches_ch = bed_ch
        .toSortedList { left, right -> left[0] <=> right[0] }
        .flatMap { rows -> rows }
        .buffer(size: params.filter_batch_size, remainder: true)
        .map { batch ->
            def columns = batch.transpose()
            // Keep one size-file name per accession, but stage shared reference files only once.
            [columns[0], columns[1], columns[2], columns[3].flatten(),
                columns[4].collect { sizes -> sizes.name }, columns[4].unique()]
        }
    // Filter BED files and create corresponding big-bed files
    FILTER_TO_BIG_BED(
        bed_batches_ch,
        file("${projectDir}/assets/peaks.as"),
        file("${projectDir}/assets/fragments.as")
    )

    // Upload once, after every filter batch completes. Do not launch an upload for an empty scan.
    upload_ch = FILTER_TO_BIG_BED.out.filtered_files
        .toList()
        .filter { batches -> !batches.isEmpty() }
        .map { batches ->
            def columns = batches.transpose()
            def accessions = columns[0].flatten().sort()
            def files = columns.drop(1).collect { column ->
                column.flatten().sort { left, right -> left.name <=> right.name }
            }
            [accessions] + files
        }
    UPLOAD_FIXED_BEDS(upload_ch, params.dry_run, params.igvf_mode)
}
