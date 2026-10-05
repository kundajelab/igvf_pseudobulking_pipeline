include { dotenv } from 'plugin/nf-dotenv'

// Write the mapping from cell name to annotation that IGVF_UPLOAD needs. It only depends on the
// metadata, so unlike the rest of the pseudobulking it is made whether or not there is RNA data.
process ANNOTATION_MAPPING {
    cpus 1
    memory '4 GB'
    conda "environments/PSEUDOBULK.yaml"
    container "${dotenv('PSEUDOBULK_IMAGE')}"
    publishDir "${params.workspace}/${params.principal_analysis.replace(",", "-")}/output",
        pattern: "${annotation_mapping}",
        mode: params.publish_mode

    input:
        path(metadata_file, arity: "1")
    output:
        path(annotation_mapping), emit: cell_name_to_annotation_mapping

    script:
    annotation_mapping = "cell_name_to_annotation_mapping.tsv"
    """
    export PYTHON_GIL=1
    pseudobulk annotation-mapping \
        --metadata-loc "${metadata_file}" \
        --output "${annotation_mapping}"
    """
}
