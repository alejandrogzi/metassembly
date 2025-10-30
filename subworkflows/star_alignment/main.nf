include { STAR_ALIGN as STAR_ALIGN_1PASS  } from './modules/nf-core/star/align/main'
include { STAR_ALIGN as STAR_ALIGN_2PASS } from './modules/nf-core/star/align/main'
include { FILTER_JUNCTIONS } from './modules/custom/junctions/main'

workflow STAR_ALIGNMENT {

    take:
    reads
    index
    gtf

    main:
    ch_versions = Channel.empty()

    ch_trimmed_reads = reads
    ch_trimmed_reads
        .multiMap { meta, reads ->
            first_pass: [meta, reads]
            second_pass: [meta, reads]
        }
        .set { ch_trimmed_split }

    ch_trimmed_split.first_pass.view { "first pass: $it" }
    index.view { "index: $it" }

    ch_star_first_pass_out = STAR_ALIGN_1PASS (
        ch_trimmed_split.first_pass,
        index.map { [ [:], it ] },
        gtf.map { [ [:], it ] },
        Channel.value([]),
        params.star_ignore_sjdbgtf,
        params.star_seq_platform ?: '',
        params.star_seq_center ?: '',
        params.star_seq_library ?: '',
        params.star_machine_type ?: '',
    )

    ch_splice_junctions = ch_star_first_pass_out.spl_junc_tab

    ch_filtered_junctions = FILTER_JUNCTIONS (
        ch_splice_junctions,
        params.star_discard_non_canonical_junctions,
        params.star_filter_junction_min_read_support,
        params.star_discard_canonical_junctions
    )

    ch_junctions_file = ch_filtered_junctions.filtered_junctions
        .map { it instanceof List ? it[1] : it }

    ch_second_pass_input = ch_trimmed_split.second_pass
        .combine(ch_junctions_file)

    ch_second_pass_input.view { "first pass: $it" }

    STAR_ALIGN_2PASS (
        ch_second_pass_input.map { meta, reads, junctions -> [meta, reads] },
        index.map { [ [:], it ] },
        gtf.map { [ [:], it ] },
        ch_second_pass_input.map { meta, reads, junctions -> junctions },
        params.star_ignore_sjdbgtf,
        params.star_seq_platform ?: '',
        params.star_seq_center ?: '',
        params.star_seq_library ?: '',
        params.star_machine_type ?: '',
    )

    ch_versions = ch_versions.mix(STAR_ALIGN_2PASS.out.versions.first())

    emit:
    bams            = STAR_ALIGN_2PASS.out.bam_sorted_aligned
    junctions       = ch_junctions_file
    log_final       = STAR_ALIGN_2PASS.out.log_final
    versions   = ch_versions   // channel: [ versions.yml ]
}
