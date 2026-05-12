include { dotenv } from 'plugin/nf-dotenv'

process COLLECT_PSEUDOBULK_QC {
    cpus 2
    memory '8 GB'
    conda "environments/CALL_PEAKS.yaml"
    container "${dotenv('CALL_PEAKS_IMAGE')}"
    publishDir "${params.workspace}/${params.principal_analysis.replace(",", "-")}/output/",
        pattern: "${pseudobulk_qc}",
        mode: params.publish_mode

    input:
        path(qc_summaries, name: "qc_summaries/*", arity: "1..*")
    output:
        path(pseudobulk_qc), emit: pseudobulk_qc_out

    script:
    pseudobulk_qc = "pseudobulk_qc.tsv.gz"
    """
    find -L qc_summaries -type f -name "*.tsv" \
    |   {
            read -r first_qc
            head -n1 "\$first_qc"
            {
                tail -n +2 "\$first_qc"
                while read -r next_qc; do
                    tail -n +2 "\$next_qc"
                done
            } \
            | sort
        } \
    | bgzip -@ ${task.cpus} -o "${pseudobulk_qc}"
    """
}
