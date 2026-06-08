/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPARE_GENOME_STAR } from '../star_index/main'
include { PREPARE_DEACON_INDEX } from '../deacon_index/main'

include { TWOBIT_TO_FA } from '../../modules/custom/ucsc/twobittofa/main'
include { CHROMSIZE } from '../../modules/custom/chromsize/main'
include { BED2GTF } from '../../modules/custom/bed2gtf/main'
include { GXF2BED } from '../../modules/custom/gxf2bed/main'

include { GUNZIP as GUNZIP_FASTA } from '../../modules/custom/gunzip/main'
include { GUNZIP as GUNZIP_GTF } from '../../modules/custom/gunzip/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREPARE_INDEXES {
    take:
        genome                              // file: /path/to/genome.{2bit/fasta}
        // --- star_index ---
        annotation                          // file: /path/to/genome.{gtf/gff/bed}
        star_index_path                     // path: /path/to/star/index/
        star_ignore_gtf_for_index           // val: boolean
        // --- deacon_index ---
        index_path                          // val: path(deacon/index)
        download_index                      // val: boolean
        make_single_index                   // val: boolean

    main:
        ch_versions = Channel.empty()

        def genome_file = file(genome, checkIfExists: true)
        def genome_path = genome_file.toString()

        ch_chrom_sizes = CHROMSIZE([[:], genome_file]).chromsize.map { it[1] }

        // INFO: if fasta is .2bit or .gz, convert or uncompress it
        if (genome_path.endsWith(".2bit")) {
            ch_fasta = TWOBIT_TO_FA([[:], genome_file]).fasta.map { it[1] }
            ch_versions = ch_versions.mix(TWOBIT_TO_FA.out.versions)
        } else if (genome_path.endsWith(".gz")) {
            ch_fasta = GUNZIP_FASTA([[:], genome_file]).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
        } else {
            ch_fasta = Channel.value(genome_file)
        }

        // INFO: preparing annotation 
        if (annotation.endsWith('.gz') || annotation.endsWith('.gtf')) {
          Channel.value(file(annotation, checkIfExists: true))
            .map { it -> [ [ id: it.baseName ], it ] }
            .set { ch_gtf }

          if (annotation.endsWith('.gz')) {
             GUNZIP_GTF(ch_gtf)

             ch_gtf = GUNZIP_GTF.out.gunzip
             ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
          }

          // INFO: converting to bed
          GXF2BED(ch_gtf)
          ch_bed = GXF2BED.out.bed

          ch_versions = ch_versions.mix(GXF2BED.out.versions)
        } else if (annotation.endsWith('.bed')) {
          // INFO: converting to gtf
          Channel.value(file(annotation, checkIfExists: true))
            .map { it -> [ [ id: it.baseName ], it ] }
            .set { ch_bed }

          BED2GTF(
            ch_bed,
            Channel.of([[], []])
          )

          ch_gtf = BED2GTF.out.gtf

          ch_versions = ch_versions.mix(BED2GTF.out.versions)
        } else {
          ch_gtf = Channel.of([:])
          ch_bed = Channel.of([[], []])
        }

        PREPARE_GENOME_STAR(
            ch_fasta,
            ch_gtf,
            star_index_path,
            star_ignore_gtf_for_index,
        )

        PREPARE_DEACON_INDEX(
            index_path,
            download_index,
            make_single_index,
            ch_fasta
        )

        ch_versions = ch_versions.mix(PREPARE_GENOME_STAR.out.versions)
        ch_versions = ch_versions.mix(PREPARE_DEACON_INDEX.out.versions)

    emit:
        genome          = ch_fasta
        star_index      = PREPARE_GENOME_STAR.out.star_index
        deacon_index    = PREPARE_DEACON_INDEX.out.deacon_index
        chrom_sizes     = ch_chrom_sizes
        annotation_gtf  = ch_gtf
        annotation_bed  = ch_bed
        versions        = ch_versions
}
