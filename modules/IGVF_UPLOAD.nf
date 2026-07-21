include { dotenv } from 'plugin/nf-dotenv'

process IGVF_UPLOAD {
    secret 'IGVF_API_KEY'
    secret 'IGVF_SECRET_KEY'
    cpus 2
    memory '4 GB'
    conda "environments/IGVF_PORTAL.yaml"
    container "${dotenv('IGVF_PORTAL_IMAGE')}"
    publishDir "${params.workspace}/${params.principal_analysis.replace(",", "-")}/output/",
        pattern: "${upload_script}",
        mode: params.publish_mode
    publishDir "${params.workspace}/${params.principal_analysis.replace(",", "-")}/output/",
        pattern: "${upload_tsvs_dir}/*.tsv",
        mode: params.publish_mode

    input:
        tuple val(principal_analysis), val(dry_run), val(igvf_mode)
        path(cell_name_to_annotation_mapping, name: "./*", arity: "1")
        path(pseudobulk_counts, name: "pseudobulks/*", arity: "1..*")
        path(pseudobulk_h5ads, name: "pseudobulks/*", arity: "1..*")
        path(fragments, name: "pseudobulks/*", arity: "1..*")
        path(peaks_bigwig, name: "pseudobulks/*", arity: "1..*")
        path(peaks_narrowpeak, name: "pseudobulks/*", arity: "1..*")
        path(peaks_pval_bigwig, name: "pseudobulks/*", arity: "1..*")
        path(pseudobulk_qc, name: "pseudobulks/*", arity: "1..*")
        path(combined_pseudobulk_qc, name: "./pseudobulk_qc.tsv.gz", arity: "1")
        path(accession_qcs, name: "analysis_accession_qc_reports/*", arity: "1..*")
    output:
        tuple path(upload_script), path("${upload_tsvs_dir}/*.tsv"), emit: upload_out

    script:
    upload_script = "upload.sh"
    upload_tsvs_dir = "upload_tsvs"
    """
    # first need to restore the pseudobulk sub-folder structure

    # make function to restore each individual pseudobulk folder by pseudobulk_id
    function restore_pseudobulk_dir {
        # move files from /pseudobulks/[pseudobulk_id].[name] to /pseudobulks/[pseudobulk_id]/[name]
        local -r pseudobulk_id=\$1
        local -r folder="pseudobulks/\$pseudobulk_id"
        mkdir -p "\$folder"
        find -L "pseudobulks" -maxdepth 1 -type f -name "\$pseudobulk_id.*" \
        | while read -r file; do
            new_name=\$(basename "\$file" | sed "s/^\$pseudobulk_id\\.//")
            mv "\$file" "\$folder/\$new_name"
        done
        if [[ -f "\$folder/sorted.tsv.gz" ]]; then
            # the fragments file needs to be renamed
            mv "\$folder/sorted.tsv.gz" "\$folder/fragments.tsv.gz"
        fi
    }

    # find pseudobulk IDs by search for QC files and restore each one
    find -L pseudobulks -type f -name "*.per_cell_qc.tsv.gz" \
    | while read -r qc_file; do
        pseudobulk_id=\$(basename "\$qc_file" .per_cell_qc.tsv.gz)
        restore_pseudobulk_dir "\$pseudobulk_id"
    done
    # ensure there are no remaining files in the pseudobulks folder
    num_remaining_files=\$(find -L pseudobulks -maxdepth 1 -type f | wc -l)
    if [[ "\$num_remaining_files" -ne 0 ]]; then
        echo "Error: found \$num_remaining_files in pseudobulks folder not assigned to any pseudobulk ID."
        exit 1
    fi

    # write the upload script
    igvf-portal gen-upload-script \
        . \
        --input-file-sets "${principal_analysis}" \
        ${dry_run ? "--dry-run" : ""} \
        --igvf-mode "${igvf_mode}"
    
    1>&2 echo "Running upload script:"
    ./upload.sh
    """
}
