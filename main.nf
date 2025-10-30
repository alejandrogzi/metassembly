#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPARE_INDEXES } from './subworkflows/prepare_indexes/main'
include { PREPROCESS_READS } from './subworkflows/preprocess_reads/main'
include { STAR_ALIGNMENT } from './subworkflows/star_alignment/main'

include { STAR_ALIGN as STAR_ALIGN_1PASS  } from './modules/nf-core/star/align/main'
include { STAR_ALIGN as STAR_ALIGN_2PASS } from './modules/nf-core/star/align/main'
include { FILTER_JUNCTIONS } from './modules/custom/junctions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow {

    main:

    ch_versions = Channel.empty()
    ch_linting_logs = Channel.empty()

    ch_start_index = Channel.empty()
    ch_deacon_index = Channel.empty()

    ch_indexes = PREPARE_INDEXES(
        params.fasta,
        params.star_sjdbgtf_path,
        params.star_ignore_sjdbgtf,
        params.star_index_path,
        params.deacon_index_path,
        params.deacon_download_index,
        params.deacon_make_single_index,
        params.deacon_multi_index_additional_genome_paths
    )

    ch_star_index = ch_indexes.star_index
    ch_star_sjdbgtf = ch_indexes.star_sjdbgtf
    ch_deacon_index = ch_indexes.deacon_index

    ch_fastqs = Channel
    .fromFilePairs("${params.input_dir}/*{1,2}.f*q.gz", checkIfExists: true, size: -1)
    .map { id, reads ->
        [
            [
                id: id,
                single_end: reads.size() == 1
            ],
            reads
        ]
    }

    ch_processed_reads = PREPROCESS_READS(
       ch_fastqs,
       ch_deacon_index
    )

    ch_processed_reads.processed_reads.view()

    ch_alignment = STAR_ALIGNMENT(
        ch_processed_reads.processed_reads,
        ch_star_index,
        ch_star_sjdbgtf,
    )
}
