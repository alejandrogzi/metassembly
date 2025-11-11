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
include { MULTIQC } from './modules/nf-core/multiqc/main'
include { EMAIL_RESULTS } from './modules/custom/email/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow METASSEMBLE {
    main:
        ch_versions = Channel.empty()
        ch_linting_logs = Channel.empty()
        ch_multiqc_files = Channel.empty()

        ch_start_index = Channel.empty()
        ch_deacon_index = Channel.empty()

        ch_indexes = PREPARE_INDEXES(
            params.fasta,
            params.star_gtf_path,
            params.star_index_path,
            params.star_ignore_gtf_for_index,
            params.deacon_index_path,
            params.deacon_download_index,
            params.deacon_make_single_index,
            params.deacon_multi_index_additional_genome_paths,
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

        ch_processed_keys = ch_processed_reads.processed_reads.map { meta, _ -> [meta, true] }
        ch_fastqs_kept = ch_fastqs
                .join(ch_processed_keys)
                .map { meta, reads, _ -> [meta, reads] }

        ch_final_reads = ch_processed_reads.processed_reads
            .combine(ch_indexes.star_index)
            .map { meta, reads, index ->
                [meta, reads, index]
            }

        ch_multiqc_files = ch_multiqc_files.mix(ch_processed_reads.fastp_json.map { it[1] })

        ch_alignment = STAR_ALIGNMENT(
            ch_final_reads,
            ch_indexes.star_gtf
        )

        ch_multiqc_files = ch_multiqc_files.mix(ch_alignment.log_final.map { it[1] })

        if (!params.skip_multiqc) {
            ch_multiqc_config        = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()
            ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()
            ch_multiqc_logo          = params.multiqc_logo   ? Channel.fromPath(params.multiqc_logo)   : Channel.empty()

            MULTIQC (
                ch_multiqc_files.collect(),
                ch_multiqc_config.toList(),
                ch_multiqc_custom_config.toList(),
                ch_multiqc_logo.collect().toList(),
                [],
                [],
            )

            ch_multiqc_report = MULTIQC.out.report
            ch_versions = ch_versions.mix(MULTIQC.out.versions)
        }

        ch_beaver = ASSEMBLY(
            ch_alignment.bams
        )

        ch_fastqs_kept
          .join(ch_alignment.bam_size, failOnMismatch: true)                        // (meta, reads, bam_size)
          .join(ch_alignment.percent_mapped, failOnMismatch: true)                  // (+ pct)
          .join(ch_processed_reads.deacon_discarded_seqs, failOnMismatch: true)     // (+ kept)
          .join(ch_processed_reads.num_trimmed_reads, failOnMismatch: true)         // (+ num_trimmed_reads)
          .join(ch_processed_reads.num_trimmed_reads_percent, failOnMismatch: true) // (+ num_trimmed_reads_percent)
          .join(ch_beaver.counts, failOnMismatch: true)                             // (+ assembled_count)
          .map {
                meta,
                reads,
                bam_size_bytes,
                pct,
                kept,
                reads_after_trim,
                reads_after_trim_percent,
                assembled_count
                ->
              def fastq_1 = file(reads[0].toUriString()).baseName
              def fastq_2 = reads.size() > 1 ? file(reads[1].toUriString()).baseName : ''

              def bam_size = (bam_size_bytes ?: 0) as long

              "${meta.id},${fastq_1},${fastq_2},${reads_after_trim},${reads_after_trim_percent},${kept ?: ''},${pct ?: ''},${bam_size / 100000000},${assembled_count ?: ''}"
          }
          .collectFile(
            name: 'samplesheet.csv',
            storeDir: "${params.outdir}/samplesheets",
            newLine: true,
          )
          .set { ch_samplesheet }

    emit:
        fastqs = ch_fastqs
        bams = ch_alignment.bams
        junctions = ch_alignment.junctions
        percent_mapped = ch_alignment.percent_mapped
        samplesheet = ch_samplesheet
        multiqc_report = ch_multiqc_report
        versions = ch_alignment.versions
}

workflow PIPELINE_COMPLETION {

    take:
    email
    email_on_fail
    plaintext_email
    outdir
    use_mailx
    ch_samplesheet

    main:

    if (params.sent_email) {
        EMAIL_RESULTS (
            email,
            email_on_fail,
            plaintext_email,
            outdir,
            use_mailx,
            ch_samplesheet
        )
    }

    workflow.onError {
        log.error "ERROR: Pipeline failed. Please refer to github issues: https://github.com/alejandrogzi/metaassembler/issues"
    }

    workflow.onComplete {
        log.info "\nPipeline completed successfully!"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    METASSEMBLE ()

    PIPELINE_COMPLETION (
        params.email_to,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.use_mailx,
        METASSEMBLE.out.samplesheet
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
