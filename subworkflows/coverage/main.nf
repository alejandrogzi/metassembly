/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BEDGRAPHTOBIGWIG } from '../../modules/custom/bigtools/bedgraphtobigwig/main'
include { WIGGLETOOLS_MEDIAN as WIGGLETOOLS } from '../../modules/custom/wiggletools/median/main'
include { UCSC_WIGTOBIGWIG as WIGTOBIGWIG } from '../../modules/nf-core/ucsc/wigtobigwig/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow COVERAGE {
    take:
        bedgraph                             // channel: [ val(meta), bedgraph ]
        chrom_sizes                          // channel: [ chrom_sizes ]

    main:
        ch_versions = Channel.empty()

        ch_bigwig = BEDGRAPHTOBIGWIG(bedgraph, chrom_sizes).bigwig

        // Collect all junction files into a single item
        ch_all_bigwig = ch_bigwig
            .map { meta, bigwig -> bigwig } // Extract just the files
            .collect() // Collect all files into a single list
            .map { files -> [ [id: 'merged_bigwigs'], files ] } // Add a generic meta map

        ch_wig = WIGGLETOOLS(ch_all_bigwig).wig
        ch_joined_bw = WIGTOBIGWIG(ch_wig, chrom_sizes).bigwig

        ch_versions = ch_versions.mix(BEDGRAPHTOBIGWIG.out.versions)
        ch_versions = ch_versions.mix(WIGGLETOOLS.out.versions)
        ch_versions = ch_versions.mix(WIGTOBIGWIG.out.versions)

    emit:
        bigwig = ch_joined_bw
        versions = ch_versions
}
