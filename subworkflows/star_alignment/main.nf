/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { STAR_ALIGN as STAR_ALIGN_1PASS } from '../../modules/nf-core/star/align/main'
include { STAR_ALIGN as STAR_ALIGN_2PASS } from '../../modules/nf-core/star/align/main'
include { JOIN_JUNCTIONS } from '../../modules/custom/junctions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow STAR_ALIGNMENT {
    take:
        reads
        gtf

    main:
        ch_versions = Channel.empty()

        ch_trimmed_reads = reads
        ch_trimmed_reads
            .multiMap { meta, reads, index ->
                first_pass: [meta, reads, index]
                second_pass: [meta, reads, index]
            }
            .set { ch_trimmed_split }

        gtf
            .multiMap { g ->
                first_pass: g
                second_pass: g
            }
            .set { gtf_split }

        if (params.star_ignore_gtf_for_mapping) {
            ch_star_first_pass_out = STAR_ALIGN_1PASS(
                ch_trimmed_split.first_pass,
                Channel.value([]),
                Channel.value([]),
                params.star_ignore_sjdbgtf,
                params.star_seq_platform ?: '',
                params.star_seq_center ?: '',
                params.star_seq_library ?: '',
                params.star_machine_type ?: ''
            )
        } else {
            ch_star_first_pass_out = STAR_ALIGN_1PASS(
                ch_trimmed_split.first_pass,
                gtf_split.first_pass,
                Channel.value([]),
                params.star_ignore_sjdbgtf,
                params.star_seq_platform ?: '',
                params.star_seq_center ?: '',
                params.star_seq_library ?: '',
                params.star_machine_type ?: ''
            )
        }

        ch_splice_junctions = ch_star_first_pass_out.spl_junc_tab

        // Collect all junction files into a single item
        ch_all_junctions = ch_splice_junctions
            .map { meta, junctions -> junctions } // Extract just the files
            .collect() // Collect all files into a single list
            .map { files -> [ [id: 'merged_junctions'], files ] } // Add a generic meta map

        ch_filtered_junctions = JOIN_JUNCTIONS(
            ch_all_junctions,
            params.junction_min_junction_length,
            params.junction_min_read_coverage
        )

        // Now ch_filtered_junctions contains a single merged junction file
        ch_junctions_file = ch_filtered_junctions.filtered_junctions
            .map { meta, junctions -> junctions } // Extract just the file

        ch_second_pass_input = ch_trimmed_split.second_pass
            .combine(ch_junctions_file)

        if (params.star_ignore_gtf_for_mapping) {
            ch_star_second_pass_out = STAR_ALIGN_2PASS(
                ch_second_pass_input.map { meta, reads, index, _junctions -> [meta, reads, index] },
                Channel.value([]),
                ch_second_pass_input.map { _meta, _reads, _index, junctions -> junctions },
                params.star_ignore_gtf_for_mapping,
                params.star_seq_platform ?: '',
                params.star_seq_center ?: '',
                params.star_seq_library ?: '',
                params.star_machine_type ?: ''
            )
        } else {
            ch_star_second_pass_out = STAR_ALIGN_2PASS(
                ch_second_pass_input.map { meta, reads, index, _junctions -> [meta, reads, index] },
                gtf_split.second_pass,
                ch_second_pass_input.map { _meta, _reads, _index, junctions -> junctions },
                params.star_ignore_gtf_for_mapping,
                params.star_seq_platform ?: '',
                params.star_seq_center ?: '',
                params.star_seq_library ?: '',
                params.star_machine_type ?: ''
            )
        }

        ch_versions = ch_versions.mix(STAR_ALIGN_2PASS.out.versions.first())
        ch_log_final = ch_star_second_pass_out.log_final

        ch_star_second_pass_out.bam_sorted_aligned
            .join(ch_star_second_pass_out.bai)
            .set { ch_bam_sorted }

    emit:
        bams = ch_bam_sorted
        junctions = ch_junctions_file
        percent_mapped = ch_log_final.map { meta, log -> [ meta, getStarPercentMapped(params, log) ] }
        versions = ch_versions // channel: [ versions.yml ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Function that parses and returns the alignment rate from the STAR log output
//
def getStarPercentMapped(params, align_log) {
    def percent_aligned = 0
    def pattern = /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/
    align_log.eachLine { line ->
        def matcher = line =~ pattern
        if (matcher) {
            percent_aligned = matcher[0][1].toFloat()
        }
    }

    return percent_aligned
}
