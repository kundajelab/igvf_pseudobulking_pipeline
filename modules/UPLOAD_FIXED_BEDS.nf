include { dotenv } from 'plugin/nf-dotenv'

process UPLOAD_FIXED_BEDS {
    secret 'IGVF_API_KEY'
    secret 'IGVF_SECRET_KEY'
    cpus 2
    memory '4 GB'
    conda "environments/IGVF_PORTAL.yaml"
    container "${dotenv('IGVF_PORTAL_IMAGE')}"
    // Keep real uploads on a queue that does not preempt.
    queue { dry_run ? params.slurm_queue : params.non_preemptable_queue }

    input:
        tuple val(bed_accessions),
            path(bed_files, arity: '1..*'),
            path(big_bed_files, arity: '1..*'),
            path(bed_jsons, name: 'input_jsons/*', arity: '1..*'),
            path(big_bed_jsons, name: 'input_jsons/*', arity: '1..*')
        val(dry_run)
        val(igvf_mode)
    output:
        path('tabular_files.json'), emit: tabular_files

    script:
    tabular_files_json = "tabular_files.json"
    """
    cat input_jsons/* > "${tabular_files_json}"

    igvf-portal register \
        ${dry_run ? "--dry-run" : "--no-dry-run"} \
        --drop-extra-fields \
        --expect-patch \
        --igvf-mode "${igvf_mode}" \
        --profile-id tabular_file \
        --infile "${tabular_files_json}"
    """
}
