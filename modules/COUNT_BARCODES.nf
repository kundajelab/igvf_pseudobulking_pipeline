include { dotenv } from 'plugin/nf-dotenv'

process COUNT_UNIQUE_BARCODES {
    cpus 2
    memory '2 GB'
    conda "environments/CALL_PEAKS.yaml"
    container "${dotenv('CALL_PEAKS_IMAGE')}"

    input:
        path(fragments, arity: "1")
        path(metadata_file, arity: "1")
    output:
        tuple path(fragments), stdout, emit: fragments_and_counts

    script:
    // The metadata file may or may not be compressed, and bgzip exits non-zero on plain input
    // ("not a compressed file -- ignored") rather than passing it through.
    metadata_reader = metadata_file.name.endsWith('.gz')
        ? "bgzip -cd -@ ${task.cpus} \"${metadata_file}\""
        : "cat \"${metadata_file}\""
    // NOTE: FS must be set to tab. Metadata columns such as cell_description hold values containing
    // spaces, and with the default separator those would be split into several fields, shifting
    // every column after them.
    // NOTE: num_codes and num_fragments are initialised so that a fragments file with no records
    // prints "0\t0" rather than empty fields, which could not be converted to integers.
    """
    awk \
        '
        BEGIN {FS="\\t"; num_codes=0; num_fragments=0}
        # the first file is the metadata: note which barcodes are kept for pseudobulking
        NR == FNR {
            if(FNR == 1) {
                for(i=1; i<=NF; ++i) {
                    if(\$i == "barcode_sample") {barcode_idx=i}
                }
                next
            }
            keep_barcodes[\$barcode_idx]=""
            next
        }
        # the second file is the fragments: count distinct barcodes, and the fragments that
        # split-fragments will hold in memory (only those whose barcode is kept)
        {
            if(!(\$4 in codes)) {codes[\$4]=""; ++num_codes}
            if(\$4 in keep_barcodes) {++num_fragments}
        }
        END {
            if(num_fragments == 0) {
                print "Found 0 fragments to pseudobulk" | "cat >&2"
                exit 1
            }
            print num_codes "\\t" num_fragments
        }
        ' \
        <(${metadata_reader}) \
        <(bgzip -cd -@ ${task.cpus} "${fragments}")
    """
}

// Count the unique barcodes in each fragments file, and how many of its fragments will be held in
// memory, emitting tuples of (fragments_file, num_unique_barcodes, num_kept_fragments).
// NOTE: the counts have to be captured with the "stdout" qualifier, which yields the whole output as
// one String, so the process is wrapped to split that single tab-separated line and emit the two
// counts as Integers.
workflow COUNT_BARCODES {
    take:
        fragments_files
        metadata_file
    main:
        COUNT_UNIQUE_BARCODES(fragments_files, metadata_file)
    emit:
        fragments_and_counts = COUNT_UNIQUE_BARCODES.out.fragments_and_counts
            .map { fragments, counts ->
                def (num_unique_barcodes, num_kept_fragments) = counts
                    .trim()
                    .tokenize('\t')
                    .collect { count -> count as Integer }
                [fragments, num_unique_barcodes, num_kept_fragments]
            }
}
