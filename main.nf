include {
    validateParameters ;
    paramsSummaryLog ;
} from 'plugin/nf-schema'

include { DOWNLOAD_ACCESSION_FILES } from './modules/DOWNLOAD_ACCESSION_FILES.nf'
include { COUNT_BARCODES } from './modules/COUNT_BARCODES.nf'
include { DETERMINE_SPECIES } from './modules/DETERMINE_SPECIES.nf'
include { SPLIT_FRAGMENTS } from './modules/SPLIT_FRAGMENTS.nf'
include { SORT_FRAGMENTS as SORT_PSEUDOREP_1 } from './modules/SORT_FRAGMENTS.nf'
include { SORT_FRAGMENTS as SORT_PSEUDOREP_2 } from './modules/SORT_FRAGMENTS.nf'
include { SORT_FRAGMENTS as SORT_PSEUDOREP_T } from './modules/SORT_FRAGMENTS.nf'
include { SORT_FRAGMENTS as SORT_SEPARATED_FRAGMENTS } from './modules/SORT_FRAGMENTS.nf'
include { MACS3 as MACS3_REP_1 } from './modules/MACS3.nf'
include { MACS3 as MACS3_REP_2 } from './modules/MACS3.nf'
include { MACS3 as MACS3_REP_T } from './modules/MACS3.nf'
include { CALL_PEAKS } from './modules/CALL_PEAKS.nf'
include { PSEUDOBULK_RNA } from './modules/PSEUDOBULK_RNA.nf'
include { SUMMARIZE_PSEUDOBULK_QC } from './modules/SUMMARIZE_PSEUDOBULK_QC.nf'
include { COMBINE_ACCESSION_QC } from './modules/COMBINE_ACCESSION_QC.nf'
include { COLLECT_PSEUDOBULK_QC } from './modules/COLLECT_PSEUDOBULK_QC.nf'
include { WRITE_SUMMARY } from './modules/WRITE_SUMMARY.nf'
include { IGVF_UPLOAD } from './modules/IGVF_UPLOAD.nf'


// Flatten channel of files and reform to channel of tuples (pseudobulk_id, file)
def flattenWithId(files_channel) {
    return files_channel
        .flatten()
        .map { file_path -> [file_path.getSimpleName(), file_path] }
}

def groupById(files_channel) {
    def sorted_list = flattenWithId(files_channel)
        .groupTuple(by: 0)
        .map { id, files -> [id, files.sort()] }
    return sorted_list
}

// Get column names present in the metadata file.
def metadataColumns(metadata_path) {
    // note: "limit: 1" stops after reading the header line.
    return file(metadata_path)
        .splitCsv(sep: '\t', limit: 1)
        .first()
        .collect { column -> column.trim() }
}

// The accessions to fetch files for, as a value channel holding a unique sorted list.
// matrix_file_accession and fragments_file_accession name individual files, so prefer them when the
// metadata has them: use both sets when both columns exist, or whichever single one exists.
// Fall back to analysis_set_accession when neither does, which leaves it to the portal to work out
// which matrix and fragments files belong to each analysis set.
// A file accession is paired with its analysis set as "<input>;<output>", because get-url names the
// downloaded file after the output accession. That keeps the "<analysis set>.bed.gz" and
// "<analysis set>.h5ad" names that SPLIT_FRAGMENTS and split_fragments.py read the accession back
// out of, and it uses the analysis set this metadata file assigns rather than re-querying the
// portal for the file's file_set.
def extractAccessions(metadata_path) {
    def columns = metadataColumns(metadata_path)
    def file_columns = ['matrix_file_accession', 'fragments_file_accession']
        .findAll { column -> column in columns }
    return channel.fromPath(metadata_path)
        .splitCsv(header: true, sep: '\t')
        .flatMap { row ->
            def analysis_set = row.analysis_set_accession?.trim()
            if (!file_columns) {
                return [analysis_set]
            }
            if (!analysis_set) {
                return []
            }
            // NOTE: coerce to String so the channel carries plain strings rather than GStrings,
            // which hold references to the row they were built from. The value goes on to be
            // interpolated into a shell command line.
            return file_columns
                .collect { column -> row[column]?.trim() }
                .findAll { accession -> accession }
                .collect { accession -> "${accession};${analysis_set}".toString() }
        }
        .filter { accession -> accession }
        .unique()
        .toSortedList()
}

workflow {
    // Print an unmistakable banner as the last thing in the log. submit-pipeline.sh sets
    // NXF_ANSI_LOG=false so the log stays readable in a file, which leaves it as a flat stream of
    // process notices: without this it is hard to tell a run that finished from one that died, or
    // one that is still going. status-pipeline.sh reports this block back.
    // NOTE: assigned here inside the entry workflow rather than declared in nextflow.config, where
    // workflow handlers are deprecated, or as a top-level "workflow.onComplete {}" statement, which
    // nextflow lint rejects as a statement mixed with script declarations. The "onComplete:" entry
    // workflow section would be tidier still, but it only compiles under NXF_SYNTAX_PARSER=v2, and
    // v2 cannot resolve the nf-dotenv "dotenv()" calls that every module uses to pick its container.
    workflow.onComplete = {
        def rule = '=' * 78
        def report = [
            rule,
            "PIPELINE ${workflow.success ? 'SUCCEEDED' : 'FAILED'}  (${workflow.runName})",
            "  duration: ${workflow.duration}",
            "  tasks:    ${workflow.stats.succeedCount} succeeded, ${workflow.stats.cachedCount} cached, " +
                "${workflow.stats.failedCount} failed, ${workflow.stats.ignoredCount} ignored",
            "  exit:     ${workflow.exitStatus}",
            "  output:   ${params.workspace}/${params.principal_analysis?.replace(",", "-")}/output",
            rule,
        ]
        // NOTE: deliberately no workflow.errorMessage here. process.shell sets -x, so for a failed
        // task that field holds the whole shell trace and would bury this summary. Nextflow already
        // prints the failing process, its exit status and its work dir immediately above.
        log.info(report.join('\n'))
    }

    // Validate input parameters
    validateParameters()

    // Record the resolved parameters, work directory and profile at the top of the log. When the
    // pipeline is submitted with submit-pipeline.sh the log is read long after the fact, so it needs
    // to say what it was actually run with. Only params differing from the schema defaults are shown.
    log.info(paramsSummaryLog(workflow))

    // Get accessions from input params if specified, otherwise collect from the metadata file.
    // NOTE: both branches have to be a value channel holding the whole list, not a plain List. The
    // list is consumed with channel operators below (buffer, map), and a plain List has no such
    // methods: passing one through raises "No signature of method: ArrayList.buffer()".
    input_accessions = params.accessions == null
        ? extractAccessions(params.metadata_file)
        : channel.value(params.accessions.split(',').collect { a -> a.trim() }.unique().sort())

    // Download all the raw files for 8 accessions at a time simultaneously. This limit:
    // 1) avoids exceeding the time limit on the redirected S3 paths
    // 2) allows downstream processing of the early batches while later batches are still downloading
    DOWNLOAD_ACCESSION_FILES(input_accessions.flatten().buffer(size: 8, remainder: true), params.igvf_mode)

    // Point to the metadata file
    metadata_file = file(params.metadata_file)

    // Determine the species, which selects which genome_data reference files to use. Unless it was
    // supplied, it takes a portal lookup to find out, so the answer only exists once a process has
    // run: it cannot be read out into a plain string here, and has to stay in a channel.
    // DETERMINE_SPECIES runs exactly once, because input_accessions is a value channel carrying the
    // whole accession list as a single item. It has to see all of them at once: get-species raises
    // "Found multiple species" when the accessions disagree, which is the check that stops a run
    // that mixes human and mouse data, and that check only works on the complete list.
    // NOTE: .first() turns the process output into a value channel. Without it the reference files
    // below are queue channels holding one item each, which the first process to read them would
    // consume, leaving every other process waiting forever.
    if (params.species) {
        species = channel.value(params.species)
    } else {
        // extractAccessions pairs a file accession with its analysis set as "<input>;<output>", but
        // the portal lookup wants bare accessions.
        DETERMINE_SPECIES(
            input_accessions.map { accessions ->
                accessions.collect { accession -> accession.tokenize(';').first() }.unique()
            },
            params.igvf_mode
        )
        species = DETERMINE_SPECIES.out.species.map { name -> name.trim() }
    }

    // Point to the species-dependent reference files. These are channels rather than plain paths
    // because species is only known once DETERMINE_SPECIES has run.
    genome_data_dir = species.map { name -> "${baseDir}/genome_data/${name}" }
    blacklist = genome_data_dir.map { dir -> file("${dir}/blacklist.bed") }
    chrom_sizes = genome_data_dir.map { dir -> file("${dir}/chr_sizes.tsv") }
    gene_info = genome_data_dir.map { dir -> file("${dir}/gene_info.csv") }
    tss_tsv = genome_data_dir.map { dir -> file("${dir}/tss.tsv") }

    // Point to some needed resource files
    peaks_auto_sql = file("${workflow.projectDir}/assets/peaks.as")
    fragments_auto_sql = file("${workflow.projectDir}/assets/fragments.as")

    // Count the unique barcodes in each fragments file, so that SPLIT_FRAGMENTS can request memory
    // proportional to the number of barcodes it has to keep QC data for
    COUNT_BARCODES(DOWNLOAD_ACCESSION_FILES.out.fragments_files.flatten(), metadata_file)

    // Split fragments files by pseudobulk
    SPLIT_FRAGMENTS(
        COUNT_BARCODES.out.fragments_and_counts,
        metadata_file,
        chrom_sizes,
        tss_tsv
    )

    // Aggregate rna data by pseudobulk
    PSEUDOBULK_RNA(DOWNLOAD_ACCESSION_FILES.out.counts_matrix_files.collect(sort: true), metadata_file, gene_info)

    // Concat and sort fragments files
    // NOTE: do not need to sort pseudoreps 1 AND 2 since they are only used for peak calling and
    // macs3 can handle unsorted input. However sorting is very fast and produces consistent order,
    // so we sort anyway.
    SORT_PSEUDOREP_1(groupById(SPLIT_FRAGMENTS.out.pseudorep_1), chrom_sizes, fragments_auto_sql, false)
    SORT_PSEUDOREP_2(groupById(SPLIT_FRAGMENTS.out.pseudorep_2), chrom_sizes, fragments_auto_sql, false)
    SORT_SEPARATED_FRAGMENTS(groupById(SPLIT_FRAGMENTS.out.separated_fragments), chrom_sizes, fragments_auto_sql, true)
    SORT_PSEUDOREP_T(groupById(SPLIT_FRAGMENTS.out.pseudorep_t), chrom_sizes, fragments_auto_sql, false)

    // Use MACS to call peaks for pseudoreps
    MACS3_REP_1(SORT_PSEUDOREP_1.out.sorted_fragments_tsv, chrom_sizes, "1")
    MACS3_REP_2(SORT_PSEUDOREP_2.out.sorted_fragments_tsv, chrom_sizes, "2")
    MACS3_REP_T(SORT_PSEUDOREP_T.out.sorted_fragments_tsv, chrom_sizes, "t")

    // Group data back together for calling peaks, this channel will yield a series of tuples (one
    // for each accession x pseudobulk combo) of the form:
    // (accession_id, pseudobulk_id, sorted_fragments, rep_1_peaks, rep_1_ppois, rep_2_peaks,
    //  rep_2_ppois, rep_t_peaks, rep_t_ppois)
    peaks_in_ch = SORT_SEPARATED_FRAGMENTS.out.sorted_fragments_tsv
        .join(MACS3_REP_1.out.output, by: 0)
        .join(MACS3_REP_2.out.output, by: 0)
        .join(MACS3_REP_T.out.output, by: 0)
    CALL_PEAKS(peaks_in_ch, chrom_sizes, blacklist, peaks_auto_sql)

    // create a channel with pseudobulked RNA files grouped together with metadata
    // of form (
    //   pseudobulk_id,
    //   rna_qc, pseudobulk_counts,
    //   frip_per_cell, fragments_per_cell, fragments_in_peaks_per_cell
    // )
    // Note that if RNA qc is missing for some pseudobulks, the 2nd and 3rd elements of the tuple
    // will be replace by a single null, and if the ATAC qc is missing for some pseudobulks, the
    // 4th, 5th, and 6th elements of the tuple will be replaced by a single null.
    // Thus we map those cases to tuples of correct length replacing missing files with an empty
    // list, which can stand in for missing files because nextflow paths can be individual files or
    // lists of files.
    summarize_qc_in_ch = flattenWithId(PSEUDOBULK_RNA.out.pseudobulked_rna_qc_reports)
        .join(flattenWithId(PSEUDOBULK_RNA.out.pseudobulk_counts), by: 0)
        .join(CALL_PEAKS.out.per_cell_stats, by: 0, remainder: true)
        .map { tup ->
            tup[1] == null ? [tup[0], [], [], tup[2], tup[3], tup[4]] :
            tup[3] == null ? [tup[0], tup[1], tup[2], [], [], []] :
            tup
        }

    // collect all the ATAC QC files and pass them as well.
    // NOTE: collect emits *nothing* for an empty channel, so if there were no fragments files (and
    // hence SPLIT_FRAGMENTS never ran) SUMMARIZE_PSEUDOBULK_QC would never be called. Use ifEmpty to
    // supply an empty list in that case, so the RNA-only QC is still summarized.
    atac_qc_files = SPLIT_FRAGMENTS.out.atac_qc_files.collect(sort: true, flat: true).ifEmpty([])
    SUMMARIZE_PSEUDOBULK_QC(summarize_qc_in_ch, atac_qc_files, metadata_file)
    // concatenate all the per-pseudobulk summary QCs and publish, keeping only the first header
    COLLECT_PSEUDOBULK_QC(SUMMARIZE_PSEUDOBULK_QC.out.qc_summary_out.collect())

    // create combined QC per accession
    rna_qc_files = PSEUDOBULK_RNA.out.all_cell_rna_qc_reports.collect(sort: true, flat: true).ifEmpty([])
    COMBINE_ACCESSION_QC(atac_qc_files, rna_qc_files, metadata_file)

    // NOTE: ifEmpty is needed on these counts because collect emits nothing at all when there were
    // no files of the corresponding type, which would stop WRITE_SUMMARY from ever being called.
    def frag_accessions = DOWNLOAD_ACCESSION_FILES.out.fragments_files
        .flatten()
        .map { frag -> frag.getSimpleName() }
        .unique()
        .collect(sort: true)
        .ifEmpty([])

    def matrix_accessions = DOWNLOAD_ACCESSION_FILES.out.counts_matrix_files
        .flatten()
        .map { frag -> frag.getSimpleName() }
        .unique()
        .collect(sort: true)
        .ifEmpty([])

    // collect some stats about the data: every accession that yielded a file of either type
    // NOTE: frag_accessions and matrix_accessions each carry a single list, so concat emits those two
    // lists rather than the accessions inside them, and flatten is needed to get back to individual
    // accessions before de-duplicating. ifEmpty is needed for the same reason as above, because
    // flatten leaves the channel empty when neither file type produced anything.
    def accessions = frag_accessions
        .concat(matrix_accessions)
        .flatten()
        .unique()
        .collect(sort: true)
        .ifEmpty([])

    WRITE_SUMMARY(
        metadata_file,
        accessions,
        frag_accessions,
        matrix_accessions,
        COMBINE_ACCESSION_QC.out.accession_qcs.collect(),
        SUMMARIZE_PSEUDOBULK_QC.out.pseudobulk_qc_out.collect()
    )

    // NOTE: the ATAC-seq and RNA-seq may be missing, so be sure to handle empty lists
    // sort to help with caching
    IGVF_UPLOAD(
        [params.principal_analysis, params.igvf_dry_run, params.igvf_mode],
        metadata_file,
        PSEUDOBULK_RNA.out.cell_name_to_annotation_mapping.collect(sort: true).ifEmpty([]),
        PSEUDOBULK_RNA.out.pseudobulk_counts.collect(sort: true).ifEmpty([]),
        PSEUDOBULK_RNA.out.pseudobulk_h5ads.collect(sort: true).ifEmpty([]),
        SORT_SEPARATED_FRAGMENTS.out.sorted_fragments_tsv
            .map{_pb_id, fragments -> fragments}
            .collect(sort: true)
            .ifEmpty([]),
        SORT_SEPARATED_FRAGMENTS.out.fragments_bigbed.collect(sort: true).ifEmpty([]),
        CALL_PEAKS.out.raw_insertions_bigwig.collect(sort: true).ifEmpty([]),
        CALL_PEAKS.out.filtered_overlap_calls.collect(sort: true).ifEmpty([]),
        CALL_PEAKS.out.pvalue_bigwig.collect(sort: true).ifEmpty([]),
        CALL_PEAKS.out.peaks_bigbed.collect(sort: true).ifEmpty([]),
        SUMMARIZE_PSEUDOBULK_QC.out.pseudobulk_qc_out.collect(sort: true),
        COLLECT_PSEUDOBULK_QC.out.pseudobulk_qc_out.collect(sort: true),
        COMBINE_ACCESSION_QC.out.accession_qcs.collect(sort: true)
    )
}
