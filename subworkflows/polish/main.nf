/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GTF_REMOVE_DIRT } from '../../modules/custom/gtf/clean/main'
include { GXF2BED } from '../../modules/custom/gxf2bed/main'
include { ISOTOOLS_ORPHAN } from '../../modules/custom/isotools/orphan/main'
include { ISOTOOLS_ORPHAN as ISOTOOLS_ORPHAN_DENOVO } from '../../modules/custom/isotools/orphan/main'
include { ISOTOOLS_FUSION } from '../../modules/custom/isotools/fusion/main'

include { XLOCI_INTRON as XLOCI_EXTRACT_INTRONS } from '../../modules/custom/xloci/intron/main.nf'
include { INTRONIC as IIC_PREDICT_SPLICEOSOME } from '../../modules/custom/intronic/main.nf'
include { ISOTOOLS_CLASSIFY_INTRON } from '../../modules/custom/isotools/classify/intron/main.nf'
include { ISOTOOLS_INTRON_RETENTION } from '../../modules/custom/isotools/intron/main.nf'
include { STRIP_OCCURRENCES as STRIP_RETENTIONS } from '../../modules/custom/strip/main'
include { STRIP_OCCURRENCES as STRIP_STRONG_RTS } from '../../modules/custom/strip/main'
include { STRIP_OCCURRENCES as STRIP_WEAK_RTS } from '../../modules/custom/strip/main'
include { STRIP_OCCURRENCES as STRIP_ARTIFACTS } from '../../modules/custom/strip/main'
include { PUBLISH as PUBLISH_FINAL_TRANSCRIPTS } from '../../modules/custom/publish/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow POLISH {
    take:
        metassembly               // channel: [ meta, gtf/bed ] 
        genome                    // channel: [ fasta ]
        annotation                // channel: [ meta, gtf ]
        repeats                   // path: /path/to/repeats.{bed/gff/gtf}
        splice_scores_dir         // path: /path/to/splice/scores/dir

    main:
        ch_versions = Channel.empty()
        ch_genome = genome.map { genome -> [ [id:genome.baseName], genome ] }

        ch_metassembly = metassembly.branch { meta, f ->
            gtf: f.name.endsWith('.gtf')
            bed: f.name.endsWith('.bed')
        }

        GTF_REMOVE_DIRT ( ch_metassembly.gtf )
        GXF2BED         ( GTF_REMOVE_DIRT.out.gtf )

        // only one branch carries data per run; the empty one spawns no tasks
        ch_cleaned_transcripts = GXF2BED.out.bed.mix( ch_metassembly.bed )

        ISOTOOLS_FUSION(
                ch_cleaned_transcripts,
                annotation
         )

        ch_splice_scores = splice_scores_dir ? Channel.fromPath(splice_scores_dir, checkIfExists: true)
          .map { it -> [ [ id: it.baseName ], it ] } : Channel.of([[], []])

        ch_full_length_transcripts = Channel.empty()
        if (annotation) {
            ISOTOOLS_ORPHAN(
                ISOTOOLS_FUSION.out.pass,
                annotation,
                ch_splice_scores
            )

            ch_full_length_transcripts = ISOTOOLS_ORPHAN.out.hq
        } else {
            ISOTOOLS_ORPHAN_DENOVO(
                GXF2BED.out.bed, // INFO: because fusion depends on annotation
                Channel.of([[], []]),
                ch_splice_scores
            )

            ch_full_length_transcripts = ISOTOOLS_ORPHAN_DENOVO.out.hq
        }


        XLOCI_EXTRACT_INTRONS(ch_genome, ch_full_length_transcripts)
        IIC_PREDICT_SPLICEOSOME(XLOCI_EXTRACT_INTRONS.out.tsv)
      
        if (repeats) {
            Channel.value([
                [ id: "repeats" ],
                file(repeats, checkIfExists: true)
            ]).set { ch_repeats }
        } else {
            ch_repeats = Channel.value([[:], []])
        }

        ch_full_length_transcripts
          .map { meta, bed -> tuple(meta.id, meta, bed) }
          .join(
            IIC_PREDICT_SPLICEOSOME.out.iic
              .map { meta, iic -> tuple(meta.id, meta, iic) }
          )
          .map { id, read_meta, bed, iic_meta, iic ->
            tuple(read_meta, bed, iic)
          }
          .set { ch_classify_inputs }

        ISOTOOLS_CLASSIFY_INTRON(
          ch_classify_inputs,
          ch_genome,
          annotation,
          ch_repeats,
          ch_splice_scores
        )

        ISOTOOLS_INTRON_RETENTION(
          ch_full_length_transcripts,
          ISOTOOLS_CLASSIFY_INTRON.out.tsv
        )

       STRIP_RETENTIONS(
          ch_full_length_transcripts,
          ISOTOOLS_INTRON_RETENTION.out.descriptor,
          "RETENTION"
        )
       STRIP_STRONG_RTS(
          STRIP_RETENTIONS.out.hq.map { meta, bed -> [ [ id: meta.id + '_no_retentions' ], bed ] },
          ISOTOOLS_INTRON_RETENTION.out.descriptor,
          "HAS_STRONG_RT"
       )
       STRIP_WEAK_RTS(
          STRIP_STRONG_RTS.out.hq.map { meta, bed -> [ [ id: meta.id + '_no_strong_rts' ], bed ] },
          ISOTOOLS_INTRON_RETENTION.out.descriptor,
          "HAS_WEAK_RT"
       )
       STRIP_ARTIFACTS(
          STRIP_WEAK_RTS.out.hq.map { meta, bed -> [ [ id: meta.id + '_no_weak_rts' ], bed ] },
          ISOTOOLS_INTRON_RETENTION.out.descriptor,
          "HAS_ARTIFACT"
       ) 

       PUBLISH_FINAL_TRANSCRIPTS(
          STRIP_ARTIFACTS.out.hq.map { meta, bed -> [ [ id: meta.id + '_clean' ], bed ] },
       )

    emit:
        hq             = STRIP_ARTIFACTS.out.hq
        retentions     = STRIP_RETENTIONS.out.discard
        strong_rts     = STRIP_STRONG_RTS.out.hq
        weak_rts       = STRIP_WEAK_RTS.out.hq
        artifacts      = STRIP_ARTIFACTS.out.hq
        introns        = ISOTOOLS_CLASSIFY_INTRON.out.tsv
        orphans        = ISOTOOLS_ORPHAN.out.scraps
        fusions        = ISOTOOLS_FUSION.out.fusion
        versions       = ch_versions
}
