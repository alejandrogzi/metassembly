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
include { ISOTOOLS_CLASSIFY_INTRON as ISOTOOLS_CLASSIFY_INTRON_TWOPASS } from '../../modules/custom/isotools/classify/intron/main.nf'
include { ISOTOOLS_INTRON_RETENTION } from '../../modules/custom/isotools/intron/main.nf'
include { ISOTOOLS_INTRON_RETENTION as ISOTOOLS_INTRON_RETENTION_TWOPASS } from '../../modules/custom/isotools/intron/main.nf'
include { STRIP_OCCURRENCES as STRIP_RETENTIONS } from '../../modules/custom/strip/main'
include { STRIP_OCCURRENCES as STRIP_RETENTIONS_TWOPASS } from '../../modules/custom/strip/main'
include { STRIP_OCCURRENCES as STRIP_STRONG_RTS } from '../../modules/custom/strip/main'
include { STRIP_OCCURRENCES as STRIP_WEAK_RTS } from '../../modules/custom/strip/main'
include { STRIP_OCCURRENCES as STRIP_ARTIFACTS } from '../../modules/custom/strip/main'
include { PUBLISH as PUBLISH_FINAL_TRANSCRIPTS } from '../../modules/custom/publish/main'

include { SORT_BED as SORT_BED_FL_TRANSCRIPTS } from '../../modules/custom/sort/main'
include { SORT_BED as SORT_BED_SCRAPS } from '../../modules/custom/sort/main'
include { SORT_BED as SORT_BED_FUSIONS } from '../../modules/custom/sort/main'
include { SORT_BED as SORT_BED_TWOPASS } from '../../modules/custom/sort/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow POLISH_TWOPASS {
    take:
        hq                       // channel: [ meta, bed ] clean first-pass HQ
        retentions               // channel: [ meta, bed ] first-pass retention discards
        intron_iic               // channel: [ meta, iic ] first-pass intronIC evidence
        genome                   // channel: [ meta, fasta ]
        annotation               // channel: [ meta, bed ]
        repeats                  // channel: [ meta, bed ]
        splice_scores            // channel: [ meta, bigwig directory ]

    main:
        // STRIP_OCCURRENCES appends `.discard` to the metadata id. Join on
        // the original id to reuse the IIC calculated before retention filtering.
        ch_classify_inputs = retentions
            .map { meta, bed ->
                def source_id = meta.id.replaceFirst(/\.discard$/, '')
                tuple(source_id, meta + [ id: "${source_id}_twopass" ], bed)
            }
            .join(
                intron_iic.map { meta, iic -> tuple(meta.id, meta, iic) }
            )
            .map { id, retention_meta, bed, iic_meta, iic ->
                tuple(retention_meta, bed, iic)
            }

        ISOTOOLS_CLASSIFY_INTRON_TWOPASS(
            ch_classify_inputs,
            genome,
            annotation,
            repeats,
            splice_scores
        )

        ch_retention_candidates = ch_classify_inputs.map { meta, bed, iic -> [meta, bed] }

        ISOTOOLS_INTRON_RETENTION_TWOPASS(
            ch_retention_candidates,
            ISOTOOLS_CLASSIFY_INTRON_TWOPASS.out.tsv
        )

        STRIP_RETENTIONS_TWOPASS(
            ch_retention_candidates,
            ISOTOOLS_INTRON_RETENTION_TWOPASS.out.descriptor,
            'RETENTION'
        )

        // First-pass HQ and its retention discards are disjoint. Mixing and
        // sorting also falls back naturally to HQ when there are no retentions.
        ch_merged_hq = hq
            .mix(STRIP_RETENTIONS_TWOPASS.out.hq)
            .map { meta, bed -> bed }
            .collectFile(
                name: "${params.prefix ?: 'polish'}.twopass.bed",
                newLine: false
            )
            .map { bed -> [ [ id: "${bed.baseName}_clean" ], bed ] }

        SORT_BED_TWOPASS(ch_merged_hq)

        ch_versions = ISOTOOLS_CLASSIFY_INTRON_TWOPASS.out.versions
            .mix(ISOTOOLS_INTRON_RETENTION_TWOPASS.out.versions)
            .mix(STRIP_RETENTIONS_TWOPASS.out.versions)
            .mix(SORT_BED_TWOPASS.out.versions)

    emit:
        hq          = SORT_BED_TWOPASS.out.sorted
        retentions  = STRIP_RETENTIONS_TWOPASS.out.discard
        introns     = ISOTOOLS_CLASSIFY_INTRON_TWOPASS.out.tsv
        versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow POLISH {
    take:
        metassembly               // channel: [ meta, gtf/bed ] 
        genome                    // channel: [ fasta ]
        annotation                // channel: [ meta, gtf ]
        repeats                   // path: /path/to/repeats.{bed/gff/gtf}
        splice_scores_dir         // path: /path/to/splice/scores/dir
        do_twopass_polish         // val: boolean

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
        ch_scraps = Channel.empty()
        if (annotation) {
            ISOTOOLS_ORPHAN(
                ISOTOOLS_FUSION.out.pass,
                annotation,
                ch_splice_scores
            )

            ch_full_length_transcripts = ISOTOOLS_ORPHAN.out.hq
            ch_scraps = ISOTOOLS_ORPHAN.out.scraps
        } else {
            ISOTOOLS_ORPHAN_DENOVO(
                GXF2BED.out.bed, // INFO: because fusion depends on annotation
                Channel.of([[], []]),
                ch_splice_scores
            )

            ch_full_length_transcripts = ISOTOOLS_ORPHAN_DENOVO.out.hq
            ch_scraps = ISOTOOLS_ORPHAN_DENOVO.out.scraps
        }

        // INFO: checkpoint to sort all beds
        SORT_BED_FL_TRANSCRIPTS(ch_full_length_transcripts)
        ch_full_length_transcripts = SORT_BED_FL_TRANSCRIPTS.out.sorted

        SORT_BED_SCRAPS(ch_scraps)
        SORT_BED_FUSIONS(ISOTOOLS_FUSION.out.fusion)

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

       STRIP_ARTIFACTS.out.hq
         .map { meta, bed -> [ [ id: meta.id + '_clean' ], bed ] }
         .set { ch_fl_hq_transcripts }

      ch_final_hq = ch_fl_hq_transcripts
      ch_final_retentions = STRIP_RETENTIONS.out.discard

      if (do_twopass_polish) {
         POLISH_TWOPASS(
            ch_fl_hq_transcripts,
            STRIP_RETENTIONS.out.discard,
            IIC_PREDICT_SPLICEOSOME.out.iic,
            ch_genome,
            annotation,
            ch_repeats,
            ch_splice_scores
         )

         ch_final_hq = POLISH_TWOPASS.out.hq
         ch_final_retentions = POLISH_TWOPASS.out.retentions
         ch_versions = ch_versions.mix(POLISH_TWOPASS.out.versions)
      }

      PUBLISH_FINAL_TRANSCRIPTS(ch_final_hq)

    emit:
        hq             = ch_final_hq
        retentions     = ch_final_retentions
        strong_rts     = STRIP_STRONG_RTS.out.hq
        weak_rts       = STRIP_WEAK_RTS.out.hq
        artifacts      = STRIP_ARTIFACTS.out.hq
        introns        = ISOTOOLS_CLASSIFY_INTRON.out.tsv
        scraps         = SORT_BED_SCRAPS.out.sorted
        fusions        = SORT_BED_FUSIONS.out.sorted
        versions       = ch_versions
}
