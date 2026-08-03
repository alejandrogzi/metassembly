/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPARE_INDEXES } from '../prepare_indexes/main'
include { PREPROCESS_READS } from '../preprocess_reads/main'
include { STAR_ALIGNMENT } from '../star_alignment/main'
include { COVERAGE } from '../coverage/main'
include { ASSEMBLY } from '../assembly/main'

include { GTF_REMOVE_DIRT } from '../../modules/custom/gtf/clean/main'
include { GXF2BED } from '../../modules/custom/gxf2bed/main'
include { ISOTOOLS_ORPHAN } from '../../modules/custom/isotools/orphan/main'
include { ISOTOOLS_ORPHAN as ISOTOOLS_ORPHAN_DENOVO } from '../../modules/custom/isotools/orphan/main'
include { ISOTOOLS_FUSION } from '../../modules/custom/isotools/fusion/main'

include { GENEPRED_LINT } from '../../modules/custom/genepred/lint/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow METASSEMBLE {
    take:
        input_dir                 // path: /path/to/input/dir
        genome                    // file: /path/to/genome.{2bit/fasta}
        annotation                // path: /path/to/annotation.{gtf/gff/bed}
        star_index_path           // path: /path/to/star/index/
        star_ignore_gtf_for_index // val: boolean
        deacon_index_path         // path: /path/to/deacon/index/
        deacon_download_index     // val: boolean
        deacon_make_single_index  // val: boolean
        star_make_coverage        // val: boolean
        skip_assembly             // val: boolean
        splice_scores_dir         // path: /path/to/splice/scores/dir
        output_dir                // path: /path/to/output/dir

    main:
        ch_versions = Channel.empty()
        ch_linting_logs = Channel.empty()
        ch_multiqc_files = Channel.empty()
        ch_multiqc_report = Channel.empty()

        ch_start_index = Channel.empty()
        ch_deacon_index = Channel.empty()

        GENEPRED_LINT(
          Channel.value(file(annotation, checkIfExists: true))
          .map { it -> [ [ id: it.baseName ], it ] }
        )

        ch_indexes = PREPARE_INDEXES(
            genome,
            annotation,
            star_index_path,
            star_ignore_gtf_for_index,
            deacon_index_path,
            deacon_download_index,
            deacon_make_single_index,
        )

        ch_fastqs = Channel
            .fromFilePairs("${input_dir}/*{1,2}.f*q.gz", checkIfExists: true, size: -1)
            .map { id, reads ->
                [
                    [
                        id: id,
                        single_end: reads.size() == 1,
                        strandedness: "unstranded"
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
            ch_indexes.annotation_gtf
        )

        if (star_make_coverage) {
          COVERAGE(
              ch_alignment.bedgraph,
              ch_indexes.chrom_sizes
          )
        }
 
        ch_metassembled_transcripts = Channel.empty()
        if (!skip_assembly) {
          ch_metassembly = ASSEMBLY(
              ch_alignment.bams,
              ch_indexes.annotation_gtf
          )

          ch_metassembly.gtf
            .map { meta, gtf -> [ [ id: gtf.baseName ], gtf ] }
            .set { ch_metassembled_transcripts }

          ch_fastqs_kept
            .join(ch_alignment.bam_size, failOnMismatch: true)                        // (meta, reads, bam_size)
            .join(ch_alignment.percent_mapped, failOnMismatch: true)                  // (+ pct)
            .join(ch_processed_reads.deacon_discarded_seqs, failOnMismatch: true)     // (+ kept)
            .join(ch_processed_reads.num_trimmed_reads, failOnMismatch: true)         // (+ num_trimmed_reads)
            .join(ch_processed_reads.num_trimmed_reads_percent, failOnMismatch: true) // (+ num_trimmed_reads_percent)
            .join(ch_metassembly.counts, failOnMismatch: true)                        // (+ assembled_count)
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
              storeDir: "${output_dir}/samplesheets",
              newLine: true,
            )
            .set { ch_samplesheet }

        } else {
          ch_samplesheet = Channel.empty()
        }

    emit:
        fastqs         = ch_fastqs
        bams           = ch_alignment.bams
        metassembly    = ch_metassembled_transcripts
        annotation_bed = ch_indexes.annotation_bed
        annotation_gtf = ch_indexes.annotation_gtf
        genome         = ch_indexes.genome
        junctions      = ch_alignment.junctions
        percent_mapped = ch_alignment.percent_mapped
        samplesheet    = ch_samplesheet
        chrom_sizes    = ch_indexes.chrom_sizes
        versions       = ch_versions
}
