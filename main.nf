include {
    validateParameters ;
    paramsHelp ;
    paramsSummaryLog ;
    samplesheetToList
} from 'plugin/nf-schema'

include { DOWNLOAD_ACCESSION_FILES } from './modules/DOWNLOAD_ACCESSION_FILES.nf'
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

workflow {
    // Validate input parameters
    validateParameters()

    // Get accessions from input params if specified, otherwise collect from the metadata file
    accessions_list = params.accessions == null
        ? channel.fromPath(params.metadata_file)
            .splitCsv(header: true, sep: '\t')
            .map { row -> row.analysis_set_accession.trim() }
            .unique()
            .toSortedList()
        : params.accessions.split(',').sort()

    // Download all the raw files for 8 accessions at a time simultaneously. This limit:
    // 1) avoids exceeding the time limit on the redirected S3 paths
    // 2) allows downstream processing of the early batches while later batches are still downloading
    DOWNLOAD_ACCESSION_FILES(accessions_list.flatten().buffer(size: 8, remainder: true))

    // Point to the metadata file
    metadata_file = file(params.metadata_file)

    // Point to the species-dependent reference files
    data_dir = "${baseDir}/genome_data/${params.species}"
    blacklist = file("${data_dir}/blacklist.bed")
    chrom_sizes = file("${data_dir}/chr_sizes.tsv")
    gene_info = file("${data_dir}/gene_info.csv")
    tss_tsv = file("${data_dir}/tss.tsv")

    // Split fragments files by pseudobulk
    SPLIT_FRAGMENTS(DOWNLOAD_ACCESSION_FILES.out.fragments_files.flatten(), metadata_file, chrom_sizes, tss_tsv)

    // Aggregate rna data by pseudobulk
    PSEUDOBULK_RNA(DOWNLOAD_ACCESSION_FILES.out.counts_matrix_files.collect(sort: true), metadata_file, gene_info)

    // Concat and sort fragments files
    // NOTE: do not need to sort pseudoreps 1 AND 2 since they are only used for peak calling and
    // macs3 can handle unsorted input. However sorting is very fast and produces consistent order,
    // so we sort anyway.
    SORT_PSEUDOREP_1(groupById(SPLIT_FRAGMENTS.out.pseudorep_1), chrom_sizes, false)
    SORT_PSEUDOREP_2(groupById(SPLIT_FRAGMENTS.out.pseudorep_2), chrom_sizes, false)
    SORT_SEPARATED_FRAGMENTS(groupById(SPLIT_FRAGMENTS.out.separated_fragments), chrom_sizes, true)
    SORT_PSEUDOREP_T(groupById(SPLIT_FRAGMENTS.out.pseudorep_t), chrom_sizes, false)

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
    CALL_PEAKS(peaks_in_ch, chrom_sizes, blacklist)

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

    // collect all the ATAC QC files and pass them as well
    atac_qc_files = SPLIT_FRAGMENTS.out.atac_qc_files.collect(sort: true, flat: true)
    SUMMARIZE_PSEUDOBULK_QC(summarize_qc_in_ch, atac_qc_files, metadata_file)
    // concatenate all the per-pseudobulk summary QCs and publish, keeping only the first header
    COLLECT_PSEUDOBULK_QC(SUMMARIZE_PSEUDOBULK_QC.out.qc_summary_out.collect())

    // create combined QC per accession
    rna_qc_files = PSEUDOBULK_RNA.out.all_cell_rna_qc_reports.collect(sort: true, flat: true)
    COMBINE_ACCESSION_QC(atac_qc_files, rna_qc_files, metadata_file)

    // collect some stats about the data
    def accessions = accessions_list.flatten().unique().collect(sort: true)
    def frag_accessions = DOWNLOAD_ACCESSION_FILES.out.fragments_files
        .flatten()
        .map { frag -> frag.getSimpleName() }
        .unique()
        .collect(sort: true)

    def matrix_accessions = DOWNLOAD_ACCESSION_FILES.out.counts_matrix_files
        .flatten()
        .map { frag -> frag.getSimpleName() }
        .unique()
        .collect(sort: true)

    WRITE_SUMMARY(
        metadata_file,
        accessions,
        frag_accessions,
        matrix_accessions,
        COMBINE_ACCESSION_QC.out.accession_qcs.collect(),
        SUMMARIZE_PSEUDOBULK_QC.out.pseudobulk_qc_out.collect()
    )

    IGVF_UPLOAD(
        [params.principal_analysis, params.igvf_dry_run, params.igvf_mode],
        PSEUDOBULK_RNA.out.cell_name_to_annotation_mapping.collect(),
        PSEUDOBULK_RNA.out.pseudobulk_counts.collect(),
        PSEUDOBULK_RNA.out.pseudobulk_h5ads.collect(),
        SORT_SEPARATED_FRAGMENTS.out.sorted_fragments_tsv.map{_pb_id, fragments -> fragments}.collect(),
        CALL_PEAKS.out.raw_insertions_bigwig.collect(),
        CALL_PEAKS.out.filtered_overlap_calls.collect(),
        CALL_PEAKS.out.pvalue_bigwig.collect(),
        SUMMARIZE_PSEUDOBULK_QC.out.pseudobulk_qc_out.collect(),
        COLLECT_PSEUDOBULK_QC.out.pseudobulk_qc_out.collect(),
        COMBINE_ACCESSION_QC.out.accession_qcs.collect()
    )
}
