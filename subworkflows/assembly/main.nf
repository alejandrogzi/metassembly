/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ALETSCH } from '../../modules/custom/aletsch/run/main'
include { BEAVER } from '../../modules/custom/beaver/run/main'
include { RENAME_GTF } from '../../modules/custom/rename/gtf/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ASSEMBLY {
    take:
        ch_bams // channel: [ val(meta), path(bam), path(bai) ]

    main:
        ch_versions = Channel.empty()

        ch_aletsch = ALETSCH(
            ch_bams
        )

        ch_renamed_gtf = RENAME_GTF(
            ch_aletsch.gtf
        )

        ch_aletsch_gtfs = ch_renamed_gtf.gtf
            .map { meta, gtf -> gtf }
            .collect()

        BEAVER(
            ch_aletsch_gtfs
        )

        ch_versions = ch_versions.mix(ALETSCH.out.versions)
        ch_versions = ch_versions.mix(BEAVER.out.versions)

    emit:
        gtf = BEAVER.out.gtf
        features = BEAVER.out.csv
        counts = ch_aletsch.assembled_transcripts
        versions = ch_versions // channel: [ versions.yml ]
}
