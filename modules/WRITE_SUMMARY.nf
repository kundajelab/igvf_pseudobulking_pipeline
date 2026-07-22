include { dotenv } from 'plugin/nf-dotenv'

process WRITE_SUMMARY {
    cpus 2
    memory '24 GB'
    conda "environments/VISUALIZE_QC.yaml"
    container "${dotenv('VISUALIZE_QC_IMAGE')}"
    publishDir "${params.workspace}/${params.principal_analysis.replace(",", "-")}/output",
        pattern: "${qc_report}",
        mode: params.publish_mode

    input:
        path(metadata_file)
        val(accession_ids)
        val(fragment_accession_ids)
        val(matrix_accession_ids)
        path(analysis_accession_qc, name: "analysis_accession_qc_reports/*", arity: "1..*")
        path(pseudobulk_qc, name: "pseudobulk_qc_reports/*", arity: "1..*")
    output:
        path(qc_report)

    script:
    input_summary_stats = "${metadata_file.getSimpleName()}.summary-stats.tsv"
    qc_report = "${metadata_file.getSimpleName()}.qc.html"
    num_frag = fragment_accession_ids.size()
    num_matrix = matrix_accession_ids.size()
    num_both = fragment_accession_ids.toSet().intersect(matrix_accession_ids).size()
    num_only_frag = fragment_accession_ids.toSet().minus(matrix_accession_ids).size()
    num_only_matrix = matrix_accession_ids.toSet().minus(fragment_accession_ids).size()
    num_neither = accession_ids.toSet().minus(fragment_accession_ids).minus(matrix_accession_ids).size()
    """
    # 1. Make the summary stats of input files
    {
        printf "accession\\tfrag\\tmatrix\\tonly-frag\\tonly-matrix\\tboth\\tneither\\n"
        printf "summary\\t%d\\t%d\\t%d\\t%d\\t%d\\t%d\\n" "$num_frag" "$num_matrix" "$num_only_frag" "$num_only_matrix" "$num_both" "$num_neither"
        awk '
            FNR == 1 { ++file_num }
            file_num == 1 { is_frag[\$1]="" }
            file_num == 2 { is_matrix[\$1]="" }
            file_num == 3 {
                printf "%s\\t%d\\t%d\\t%d\\t%d\\t%d\\t%d\\n", \$1, \$1 in is_frag, \$1 in is_matrix, \
                    (\$1 in is_frag) && !(\$1 in is_matrix), (\$1 in is_matrix) && !(\$1 in is_frag), \
                    (\$1 in is_frag) && (\$1 in is_matrix), !(\$1 in is_frag) && !(\$1 in is_matrix)
            }' \
            <(echo -e "${fragment_accession_ids.join('\\n')}") \
            <(echo -e "${matrix_accession_ids.join('\\n')}") \
            <(echo -e "${accession_ids.join('\\n')}")
    } > "${input_summary_stats}"


    # 2. create a temporary directory in case pandas/pola.rs are writing temporary files    
    export TMPDIR=\$(mktemp -d -p .)
    trap 'rm -rf "\$TMPDIR"' EXIT

    # 3. call visualize-qc to make the report
    export PYTHONFAULTHANDLER=1
    visualize-qc \
        --output "${qc_report}" \
        --table-qc "${input_summary_stats}" \
        --accession-qc analysis_accession_qc_reports \
        --pseudobulk-qc pseudobulk_qc_reports
    """
}
