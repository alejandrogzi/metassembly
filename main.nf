#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPARE_INDEXES } from './subworkflows/prepare_indexes/main'
include { PREPROCESS_READS } from './subworkflows/preprocess_reads/main'
include { STAR_ALIGNMENT } from './subworkflows/star_alignment/main'
include { ASSEMBLY } from './subworkflows/assembly/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow METASSEMBLE {
    main:
        ch_versions = Channel.empty()
        ch_linting_logs = Channel.empty()

        ch_start_index = Channel.empty()
        ch_deacon_index = Channel.empty()

        ch_indexes = PREPARE_INDEXES(
            params.fasta,
            params.star_gtf_path,
            params.star_ignore_sjdbgtf,
            params.star_index_path,
            params.deacon_index_path,
            params.deacon_download_index,
            params.deacon_make_single_index,
            params.deacon_multi_index_additional_genome_paths
        )

        ch_fastqs = Channel
            .fromFilePairs("${params.input_dir}/*{1,2}.f*q.gz", checkIfExists: true, size: -1)
            .map { id, reads ->
                [
                    [
                        id: id,
                        single_end: reads.size() == 1,
                        strandedness: "paired_end"
                    ],
                    reads
                ]
            }

        ch_processed_reads = PREPROCESS_READS(
            ch_fastqs,
            ch_indexes.deacon_index
        )

        ch_final_reads = ch_processed_reads.processed_reads
            .combine(ch_indexes.star_index)
            .map { meta, reads, index ->
                [meta, reads, index]
            }

        ch_alignment = STAR_ALIGNMENT(
            ch_final_reads,
            ch_indexes.star_gtf
        )

        ch_beaver = ASSEMBLY(
            ch_alignment.bams
        )

        // Create channel with original input info and BAM paths
        // ch_alignment.bams.map { meta, reads, bai -> [ meta.id, meta, reads, bai ] }
        //     .join(ch_alignment.percent_mapped)
        //     .join(ch_processed_reads.deacon_discarded_seqs)
        //     .transpose()
        //     .map { id, bam_meta, reads, meta, percent_mapped, deacon_stats ->
        //         def fastq_1 = reads[0].toUriString()
        //         def fastq_2 = reads.size() > 1 ? reads[1].toUriString() : ''
        //         def mapped = percent_mapped != null ? percent_mapped : ''
        //         def deacon_stats = deacon_stats != null ? deacon_stats : ''

        //         return "${meta.id},${fastq_1},${fastq_2},${mapped},${deacon_stats}"
        //     }
        //     .collectFile(
        //         name: 'samplesheet_with_bams.csv',
        //         storeDir: "${params.outdir}/samplesheets",
        //         newLine: true,
        //         seed: 'sample,fastq_1,fastq_2,percent_mapped'
        //     )

    emit:
        fastqs = ch_fastqs
        bams = ch_alignment.bams
        junctions = ch_alignment.junctions
        percent_mapped = ch_alignment.percent_mapped
        versions = ch_alignment.versions
}


// workflow PIPELINE_COMPLETION {

//     take:
//     email           //  string: email address
//     email_on_fail   //  string: email address sent on pipeline failure
//     plaintext_email // boolean: Send plain-text email instead of HTML
//     outdir          //    path: Path to output directory where results will be published
//     monochrome_logs // boolean: Disable ANSI colour codes in log output
//     hook_url        //  string: hook URL for notifications
//     map_status         // map: pass/fail status per sample for mapping

//     main:
//     def pass_mapped_reads  = [:]

//     summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

//     //
//     // Completion email and summary
//     //
//     workflow.onComplete {
//         if (email || email_on_fail) {
//             completionEmail(
//                 summary_params,
//                 email,
//                 email_on_fail,
//                 plaintext_email,
//                 outdir,
//                 monochrome_logs,
//                 multiqc_reports.getVal(),
//             )
//         }

//         rnaseqSummary(monochrome_logs, pass_mapped_reads, pass_trimmed_reads, pass_strand_check)

//         if (hook_url) {
//             imNotification(summary_params, hook_url)
//         }
//     }

//     workflow.onError {
//         log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
//     }
// }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    METASSEMBLE ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
