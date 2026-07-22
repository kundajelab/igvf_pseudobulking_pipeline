include { dotenv } from 'plugin/nf-dotenv'

// Use aria2c to download all accession files in parallel
process DETERMINE_SPECIES {
    secret 'IGVF_API_KEY'
    secret 'IGVF_SECRET_KEY'
    cpus 1
    memory '8 GB'
    conda "environments/IGVF_PORTAL.yaml"
    container "${dotenv('IGVF_PORTAL_IMAGE')}"
    cache "lenient"
    maxForks 1

    input:
        val(accessions)
        val(igvf_mode)
    output:
        stdout emit: species

    script:
    // NOTE: get-species splits its argument on commas, so the accessions have to be joined rather
    // than interpolated as a list: "${accessions}" would render a groovy List as "[A, B]" and the
    // brackets would end up in the first and last accession it tries to look up.
    """
    igvf-portal get-species "${accessions.join(',')}" --igvf-mode "${igvf_mode}"
    """
}
