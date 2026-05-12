include { dotenv } from 'plugin/nf-dotenv'

process CALL_PEAKS {
    cpus 1
    memory '8 GB'
    conda "environments/CALL_PEAKS.yaml"
    container "${dotenv('CALL_PEAKS_IMAGE')}"
    publishDir "${params.workspace}/${params.principal_analysis.replace(",", "-")}/output/pseudobulks/${pseudobulk_id}",
        pattern: "${raw_insertions_bigwig}",
        saveAs: { _fileName -> "raw_insertions.bw" },
        mode: params.publish_mode
    publishDir "${params.workspace}/${params.principal_analysis.replace(",", "-")}/output/pseudobulks/${pseudobulk_id}",
        pattern: "${filtered_overlap_calls}",
        saveAs: { _fileName -> "peaks.narrowPeak.gz" },
        mode: params.publish_mode
    publishDir "${params.workspace}/${params.principal_analysis.replace(",", "-")}/output/pseudobulks/${pseudobulk_id}",
        pattern: "${pvalue_bigwig}",
        saveAs: { _fileName -> "peaks_minuslog10pval.bw" },
        mode: params.publish_mode

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
    output:
        tuple val(pseudobulk_id),
            path(frip_per_cell),
            path(fragments_per_cell),
            path(fragments_in_peaks_per_cell),
            emit: per_cell_stats
        path(raw_insertions_bigwig), emit: raw_insertions_bigwig
        path(filtered_overlap_calls), emit: filtered_overlap_calls
        path(pvalue_bigwig), emit: pvalue_bigwig

    script:
    base = pseudobulk_id
    overlap_output = "${base}.peaks_overlap_unfiltered.narrowPeak"
    filtered_overlap_calls = "${base}.peaks.narrowPeak.gz"
    combined_ppois = "${base}.combined_ppois.bdg"
    combined_sorted_ppois = "${base}.combined_ppois_sorted.bdg"
    pvalue_bigwig = "${base}.peaks_minuslog10pval.bw"
    raw_insertions_ppois = "${base}.raw_insertions.bdg"
    raw_insertions_bigwig = "${base}.raw_insertions.bw"
    fragments_per_cell = "${base}.fragments_per_cell.tsv"
    fragments_in_peaks = "${base}.fragments_in_peaks.tsv"
    fragments_in_peaks_per_cell = "${base}.fragments_in_peaks_per_cell.tsv"
    frip_per_cell = "${base}.frip_per_cell.tsv"
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
    | bedtools intersect \
        -u \
        -a stdin \
        -b "${rep_2_top_peak_calls}" \
        -g "${chrom_sizes}" \
        -f "${params.min_overlap}" \
        -F "${params.min_overlap}" \
        -e -sorted \
    > "${overlap_output}"

    1>&2 echo "Filtering blacklist peaks"
    bedtools intersect -v -a "${overlap_output}" -b "${blacklist}" \
    | bgzip -@ ${task.cpus} -o "${filtered_overlap_calls}"

    1>&2 echo "Combining p-value bedgraphs"
    macs3 cmbreps \
        -m fisher \
        -i "${rep_1_clipped_ppois}" "${rep_2_clipped_ppois}" "${rep_t_clipped_ppois}" \
        -o "${combined_ppois}"

    1>&2 echo "Sorting combined ppois"
    tail -n +2 "${combined_ppois}" \
    | sort-bed.sh "${chrom_sizes}" \
    > "${combined_sorted_ppois}"

    1>&2 echo "Converting p-value bedgraphs to bigwigs"
    bedGraphToBigWig "${combined_sorted_ppois}" "${chrom_sizes}" "${pvalue_bigwig}"

    1>&2 echo "Making raw insertion bigwig"
    bedtools genomecov -i "${rep_t_top_peak_calls}" -g "${chrom_sizes}" -bg > "${raw_insertions_ppois}"
    bedGraphToBigWig "${raw_insertions_ppois}" "${chrom_sizes}" "${raw_insertions_bigwig}"

    1>&2 echo "Computing per cell FRiP"
    bgzip -cd "${separated_fragments}" \
    | cut -f4 \
    | sort \
    | uniq -c \
    | awk -v OFS='\\t' '{print \$2, \$1}' \
    > "${fragments_per_cell}"

    bedtools intersect \
        -a "${separated_fragments}" \
        -b ${filtered_overlap_calls} \
        -u > "${fragments_in_peaks}"

    cut -f4 "${fragments_in_peaks}" \
    | sort \
    | uniq -c \
    | awk -v OFS='\\t' '{print \$2, \$1}' \
    > "${fragments_in_peaks_per_cell}"

    join \
        -a 1 \
        -1 1 \
        -2 1 \
        <(sort "${fragments_per_cell}") \
        <(sort "${fragments_in_peaks_per_cell}") \
    | awk -v OFS='\\t' '{if(\$3=="") \$3=0; frip=\$3/\$2; print \$1, frip}' \
    > "${frip_per_cell}"
    """
}
