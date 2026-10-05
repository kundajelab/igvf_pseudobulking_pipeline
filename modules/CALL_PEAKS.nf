include { dotenv } from 'plugin/nf-dotenv'

process CALL_PEAKS {
    cpus 1
    memory '8 GB'
    conda "environments/CALL_PEAKS.yaml"
    container "${dotenv('CALL_PEAKS_IMAGE')}"

    input:
        tuple val(pseudobulk_id),
            path(separated_fragments),
            path(rep_1_top_peak_calls),
            path(rep_1_clipped_ppois),
            path(rep_2_top_peak_calls),
            path(rep_2_clipped_ppois),
            path(rep_t_top_peak_calls),
            path(rep_t_clipped_ppois)
        path(chrom_sizes)
        path(blacklist)
        path(peaks_auto_sql)
    output:
        tuple val(pseudobulk_id),
            path(frip_per_cell),
            path(fragments_per_cell),
            path(fragments_in_peaks_per_cell),
            emit: per_cell_stats
        path(raw_insertions_bigwig), emit: raw_insertions_bigwig
        path(filtered_overlap_calls), emit: filtered_overlap_calls
        path(filtered_overlap_bigbed), emit: peaks_bigbed
        path(pvalue_bigwig), emit: pvalue_bigwig

    script:
    base = pseudobulk_id
    rep_t_rep_1_overlap = "${base}.peaks_overlap_t_1.narrowPeak"
    overlap_output = "${base}.peaks_overlap_unfiltered.narrowPeak"
    filtered_overlap_calls = "${base}.peaks.narrowPeak.gz"
    filtered_overlap_bigbed = "${base}.peaks.narrowPeak.bb"
    combined_ppois = "${base}.combined_ppois.bdg"
    combined_sorted_ppois = "${base}.combined_ppois_sorted.bdg"
    pvalue_bigwig = "${base}.peaks_minuslog10pval.bw"
    raw_insertions_ppois = "${base}.raw_insertions.bdg"
    raw_insertions_bigwig = "${base}.raw_insertions.bw"
    fragments_per_cell = "${base}.fragments_per_cell.tsv"
    fragments_in_peaks_per_cell = "${base}.fragments_in_peaks_per_cell.tsv"
    frip_per_cell = "${base}.frip_per_cell.tsv"
    // Cap sort's buffer at half the task's memory, leaving room for the rest of the pipeline. By
    // default sort sizes it from the node's physical memory, which can exceed the job's limit.
    sort_buffer_size = "${task.memory.toMega().intdiv(2)}M"
    """
    1>&2 echo "Intersecting peaks"
    bedtools intersect \
        -u \
        -a "${rep_t_top_peak_calls}" \
        -b "${rep_1_top_peak_calls}" \
        -g "${chrom_sizes}" \
        -f "${params.min_overlap}" \
        -F "${params.min_overlap}" \
        -e -sorted \
    > "${rep_t_rep_1_overlap}"
    # NOTE: don't pipe into the second intersect: with -sorted, bedtools stops reading -a once -b is
    # exhausted, so the first intersect can be killed by SIGPIPE (exit 141) while still writing.
    bedtools intersect \
        -u \
        -a "${rep_t_rep_1_overlap}" \
        -b "${rep_2_top_peak_calls}" \
        -g "${chrom_sizes}" \
        -f "${params.min_overlap}" \
        -F "${params.min_overlap}" \
        -e -sorted \
    > "${overlap_output}"

    1>&2 echo "Filtering blacklist peaks"
    # Use bedtools to keep peaks that don't overlap the blacklist.
    # NOTES on fixes that require awk:
    # 1) below that the portal audits `score` values over 1000, because some visualization
    #    software will error in that case. Use awk to cap score.
    # 2) MACS3 can extend intervals past the end of a genomic interval, which can cause
    #    bedToBigBed and the IGVF Portal to error, so cap END at the maximum genome coordinate
    bedtools intersect -v -a "${overlap_output}" -b "${blacklist}" \
    | awk \
        -F '\\t' \
        -v OFS='\\t' \
        '
        FNR==1 { ++file_num }
        file_num == 1 { max_end[\$1]=\$2 }
        file_num == 2 {
            if(\$5 > 1000) { \$5 = 1000 }
            if(\$3 > max_end[\$1]) { \$3 = max_end[\$1] }
            print \$0
        }
        '\
        "${chrom_sizes}" \
        - \
    | bgzip -@ ${task.cpus} -o "${filtered_overlap_calls}"

    # make bigbed version of filtered peaks
    bedToBigBed \
        -type=bed6+4 \
        -as="${peaks_auto_sql}" \
        "${filtered_overlap_calls}" \
        "${chrom_sizes}" \
        "${filtered_overlap_bigbed}"

    1>&2 echo "Combining p-value bedgraphs"
    macs3 cmbreps \
        -m fisher \
        -i "${rep_1_clipped_ppois}" "${rep_2_clipped_ppois}" "${rep_t_clipped_ppois}" \
        -o "${combined_ppois}"

    1>&2 echo "Sorting combined ppois"
    tail -n +2 "${combined_ppois}" \
    | sort-bed.sh --buffer-size "${sort_buffer_size}" "${chrom_sizes}" \
    > "${combined_sorted_ppois}"

    1>&2 echo "Converting p-value bedgraphs to bigwigs"
    bedGraphToBigWig "${combined_sorted_ppois}" "${chrom_sizes}" "${pvalue_bigwig}"

    1>&2 echo "Making raw insertion bigwig"
    bedtools genomecov -i "${rep_t_top_peak_calls}" -g "${chrom_sizes}" -bg > "${raw_insertions_ppois}"
    bedGraphToBigWig "${raw_insertions_ppois}" "${chrom_sizes}" "${raw_insertions_bigwig}"

    1>&2 echo "Computing per cell FRiP"
    # Count fragments per barcode with an awk hash rather than sorting every fragment: the number of
    # barcodes is small, so only the per-barcode counts are sorted (to make the output repeatable).
    function count_per_barcode {
        awk -F '\\t' -v OFS='\\t' '
            { ++count[\$4] }
            END { for(barcode in count) { print barcode, count[barcode] } }
        ' \
        | LC_ALL=C sort -k1,1
    }

    bgzip -cd "${separated_fragments}" \
    | count_per_barcode \
    > "${fragments_per_cell}"

    bedtools intersect \
        -a "${separated_fragments}" \
        -b "${filtered_overlap_calls}" \
        -u \
    | count_per_barcode \
    > "${fragments_in_peaks_per_cell}"

    # NOTE: test FILENAME rather than FNR==NR, because the in-peaks counts may be empty.
    awk \
        -F '\\t' \
        -v OFS='\\t' \
        -v in_peaks_file="${fragments_in_peaks_per_cell}" \
        '
        FILENAME == in_peaks_file { in_peaks[\$1] = \$2; next }
        { print \$1, (\$1 in in_peaks ? in_peaks[\$1] : 0) / \$2 }
        ' \
        "${fragments_in_peaks_per_cell}" \
        "${fragments_per_cell}" \
    > "${frip_per_cell}"
    """
}
