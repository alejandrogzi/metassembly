/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { XORF } from '../../modules/xorf/src/subworkflows/xorf/main.nf'
include { WGET as WGET_SAMBA_WEIGHTS } from '../../modules/xorf/src/modules/wget/main.nf'
include { WGET as WGET_PROTEIN_DATABASE } from '../../modules/xorf/src/modules/wget/main.nf'
include { UNTAR } from '../../modules/xorf/src/modules/untar/main.nf'
include { GUNZIP as GUNZIP_DATABASE } from '../../modules/xorf/src/modules/gunzip/main.nf'
include { FASTA_MERGE } from '../../modules/xorf/src/modules/diamond/merge/main.nf'
include { DIAMOND_MAKEDB } from '../../modules/xorf/src/modules/diamond/makedb/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow XORF_RUN {
    take:
        regions                 // channel: [ meta, bed ] first-pass HQ transcripts
        sequence                // channel: [ meta, fasta ] genome (2bit ok)

    main:
        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            SAMBA WEIGHTS
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        ch_samba_weights = Channel.empty()
        if (params.xorf_samba_local_weights) {
            ch_samba_weights = Channel.value(
              file(params.xorf_samba_local_weights, checkIfExists: true)
            ).map { path -> [ [id : path.baseName ], path ] }
        } else {
            WGET_SAMBA_WEIGHTS(
              Channel.value(
                params.xorf_samba_weights
              ).map { url -> [ [id : url.tokenize('/')[-1]], url ] }
            )
            ch_samba_weights = WGET_SAMBA_WEIGHTS.out.outfile
        }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            PROTEIN DATABASE
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        // Mirrors the standalone xorf database preparation (xorf src/main.nf).
        // BLAST declares `each path(database)`, so ch_database carries bare files.
        ch_database = Channel.empty()
        ch_database_versions = Channel.empty()
        if (params.xorf_custom_database) {
            if (params.xorf_custom_database.endsWith('.dmnd')) {
                ch_database = Channel.fromPath(params.xorf_custom_database, checkIfExists: true)
            } else if (params.xorf_custom_database.endsWith('.dmnd.gz')) {
                GUNZIP_DATABASE(
                    Channel.value(
                        [ [id: params.xorf_custom_database.tokenize('/')[-1]], params.xorf_custom_database ]
                    )
                )
                GUNZIP_DATABASE.out.gunzip
                    .map { meta, it -> it }
                    .set { ch_database }
            } else if (params.xorf_custom_database.endsWith('.fa') || params.xorf_custom_database.endsWith('.fasta')
                        || params.xorf_custom_database.endsWith('.fa.gz') || params.xorf_custom_database.endsWith('.fasta.gz')) {
                WGET_PROTEIN_DATABASE(
                    Channel.value('https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz')
                    .map { it -> [ [ id: 'uniprot_sprot.fasta.gz' ], it ] }
                )
                FASTA_MERGE(
                    Channel.fromPath(params.xorf_custom_database, checkIfExists: true)
                        .map { it -> [ [ id: it.baseName ], it ] }
                        .combine(WGET_PROTEIN_DATABASE.out.outfile.map { meta, raw -> raw })
                )
                DIAMOND_MAKEDB(
                    FASTA_MERGE.out.fasta.map { meta, fasta -> [ [ id: 'merged_database' ], fasta ] },
                    [],
                    [],
                    []
                )
                DIAMOND_MAKEDB.out.db
                    .map { meta, it -> it }
                    .set { ch_database }
                ch_database_versions = FASTA_MERGE.out.versions.mix(DIAMOND_MAKEDB.out.versions)
            } else {
                error """
                ERROR: custom_database extension not recognized.
                Please provide a custom database in one of these formats:
                  .dmnd/.dmnd.gz   -> replaces the default database entirely
                  .fa/.fasta/.fa.gz/.fasta.gz -> appended to the default SwissProt database and reindexed
                """.stripIndent()
            }
        } else {
            WGET_PROTEIN_DATABASE(
                Channel.value('https://zenodo.org/records/21399231/files/diamond_db.tar.gz?download=1')
                    .map { it -> [ [ id: 'uniprot_sprot.tar.gz' ], it ] }
            )
            UNTAR(WGET_PROTEIN_DATABASE.out.outfile)
            ch_database = UNTAR.out.contents.map { meta, it -> it }
        }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            MAIN WORKFLOW
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        // meta.id carries the source transcript id (`<t>_clean`); XORF keeps it as the
        // group key and prefixes the renamed ORF ids with it. `chr` is only a grouping
        // hash: CHUNKER uses it as its prefix and the subworkflow groups by
        // `${meta.id}@${meta.chr}`. Using meta.id keeps every XORF-derived filename
        // deterministic (required by the e2e golden checks).
        XORF(
            regions.map { meta, f -> [ meta + [ chr: meta.id ], f ] },
            sequence.map { meta, fa -> fa },   // XORF takes a bare-file sequence channel
            ch_database,
            "${params.output_dir}/09_polish/xorf",
            params.xorf_chunk_size ?: 20,
            ch_samba_weights,
            params.xorf_predict_keep_raw,
            params.xorf_selenocysteine_sites,                          
            params.xorf_skip_netstart,
            params.xorf_rename_deactivate,
            params.xorf_do_polishing,
            params.xorf_skip_joined_concat,
            false, null, null,             // run_only_on / run_only_mode / run_only_target — masking not ported
            ch_database_versions
        )

    emit:
        files    = XORF.out.files     // [ meta, bed, tsv ] renamed/merged ORF predictions
        counts   = XORF.out.counts
        versions = XORF.out.versions
}
