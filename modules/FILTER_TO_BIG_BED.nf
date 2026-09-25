include { dotenv } from 'plugin/nf-dotenv'

process FILTER_TO_BIG_BED {
    cpus 2
    memory { 10.B * bed_files.collect { bed -> bed.size() }.max() }
    conda "${moduleDir}/../environments/CALL_PEAKS.yaml"
    container "${dotenv('CALL_PEAKS_IMAGE')}"

    input:
        tuple val(bed_accessions),
            path(bed_files, arity: "1..*"),
            path(bed_jsons, arity: "1..*"),
            path(big_bed_jsons, arity: "1..*"),
            val(chr_sizes_names),
            path(chr_sizes_files, arity: "1..*")
        path(peaks_auto_sql)
        path(fragments_auto_sql)
    output:
        tuple val(bed_accessions),
            path("${output_dir}/*.bed.gz", arity: '0..*'),
            path(big_beds, arity: '1..*'),
            path("${output_dir}/*.bed.json", arity: '0..*'),
            path(fixed_big_bed_jsons, arity: '1..*'),
            emit: filtered_files

    script:
    output_dir = "output"
    big_beds = bed_accessions.collect { accession -> "${output_dir}/${accession}.bb" }
    fixed_big_bed_jsons = big_bed_jsons.collect { record -> "${output_dir}/${record.name}" }
    """
    accessions=( "${bed_accessions.join('" "')}" )
    chr_sizes_names=( "${chr_sizes_names.join('" "')}" )

    mkdir -p "${output_dir}"

    function update_record {
        local -r fixed_file="\$1"
        local -r record_json="\$2"
        local -r checksum=\$(md5sum "\$fixed_file" | awk '{print \$1}')
        # The upload process stages the data files by basename in its working directory.
        jq \
            --arg md5sum "\$checksum" \
            --arg submitted_file_name "\${fixed_file##*/}" \
            '.md5sum = \$md5sum | .submitted_file_name = \$submitted_file_name' \
            "\$record_json" > "${output_dir}/\$record_json"
    }

    for index in "\${!accessions[@]}"; do
        accession="\${accessions[\$index]}"
        chr_sizes="\${chr_sizes_names[\$index]}"
        bed_json="\${accession}.bed.json"
        content_type=\$(jq -er '.content_type' "\$bed_json")
        case "\$content_type" in
            peaks)
                bed_file="\${accession}.tsv.gz"
                bed_type=bed6+4
                auto_sql="${peaks_auto_sql}"
                ;;
            fragments)
                bed_file="\${accession}.bed.gz"
                bed_type=bed3+2
                auto_sql="${fragments_auto_sql}"
                ;;
            *)
                echo "Unsupported content_type '\$content_type' in \$bed_json" >&2
                exit 1
                ;;
        esac
        if [[ "\$(jq -r '.upload_status' "\$bed_json")" == "validated" ]]; then
            # this BED is okay, we just needed it to make the big-bed
            filtered_bed="\$bed_file"
        else
            # the bed needs to be filtered
            filtered_bed="${output_dir}/\${accession}.bed.gz"
            1>&2 echo "Filtering \$content_type \$bed_file with \$chr_sizes"
            # Apply the CALL_PEAKS interval-end filter; fragment counts are not BED scores.
            awk \
                -F '\\t' \
                -v OFS='\\t' \
                -v content_type="\$content_type" \
                '
                FNR==1 { ++file_num }
                file_num == 1 { max_end[\$1]=\$2 }
                file_num == 2 {
                    if(content_type == "peaks" && \$5 > 1000) { \$5 = 1000 }
                    if(\$3 > max_end[\$1]) { \$3 = max_end[\$1] }
                    print \$0
                }
                ' \
                "\$chr_sizes" \
                <(bgzip -cd -@ ${task.cpus} "\$bed_file") \
            | sort-bed.sh "\$chr_sizes" - \
            | bgzip -@ ${task.cpus} -o "\$filtered_bed"

            # update the bed record for later upload
            update_record "\$filtered_bed" "\$bed_json"
        fi

        1>&2 echo "Converting \$content_type \$bed_file to big-bed"
        big_bed="${output_dir}/\${accession}.bb"
        bedToBigBed \
            -type="\$bed_type" \
            -as="\$auto_sql" \
            "\$filtered_bed" \
            "\$chr_sizes" \
            "\$big_bed"

        big_bed_json="\${accession}.bb.json"
        update_record "\$big_bed" "\$big_bed_json"
    done
    """
}
