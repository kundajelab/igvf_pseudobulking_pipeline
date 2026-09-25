include { dotenv } from 'plugin/nf-dotenv'

process FIND_BAD_BEDS {
    secret 'IGVF_API_KEY'
    secret 'IGVF_SECRET_KEY'
    cpus 1
    memory '4 GB'
    container "${dotenv('IGVF_PORTAL_IMAGE')}"

    input:
        val(igvf_mode)
        val(max_pseudobulks)
        val(lab)
    output:
        path("*.bed.json"), optional: true, emit: bed_records
        path("*.bb.json"), optional: true, emit: big_bed_records
        path("*.ref.txt"), optional: true, emit: ref_accession_files

    script:
    max_pseudobulks_arg = max_pseudobulks == "" ? "" : "--max-pseudobulks ${max_pseudobulks}"
    """
    igvf-portal find-bad-beds \
        --output . \
        --igvf-mode "${igvf_mode}" \
        --lab "${lab}" \
        ${max_pseudobulks_arg}
    """
}
