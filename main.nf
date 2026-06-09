#!/usr/bin/env nextflow

/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    xasm

    meta-assemble bulk RNA-seq datasets at scale
    Authors: Alejandro Gonzales-Irribarren, Michael Hiller

    GitHub:  https://github.com/hillerlab/xasm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.help) {
    log.info """
    xasm v${workflow.manifest.version}
    Meta-assemble bulk RNA-seq datasets at scale

    Authors: ${workflow.manifest.author}
    Github:  ${workflow.manifest.homePage}

    Usage (full run):
        nextflow run main.nf \\
            --input_dir   PATH    Path to input directory (fastqs)
            --genome      PATH    Path to genome file (fasta/2bit)
            --annotation  PATH    Path to annotation file (gtf/gff/bed) [defines STAR + iso-orphan inputs]
            --splice_scores_dir  PATH      Path to splicing scores dir
            --repeats     PATH      Path to repeats.{bed/gff/gtf}
            --output_dir  PATH    Output directory [default: ./results]
            -profile         apptainer,slurm

    Pass all parameters from a JSON file (replaces old --params_from_file):
        nextflow run main.nf -params-file my_params.json

    Required parameters (full run + fill/clean):
        --input_dir           PATH    Path to input directory (fastqs)
        --genome              PATH    Path to genome file (fasta/2bit)
        --annotation          PATH    Path to annotation file (gtf/gff/bed) [defines STAR + iso-orphan inputs]
        --output_dir          PATH    Output directory [default: ./results]

    Optional parameters (common):
        --splice_scores_dir   PATH      Path to splicing scores dir 
        --repeats             PATH      Path to repeats.{bed/gff/gtf}
        --from                STRING    Checkpoint to resume from [options: polish]
        --polish_path         PATH      Path to polished assembly [required if --from polish]

    Profiles:
        local       Run on local machine (default)
        slurm       Submit jobs to SLURM cluster
        conda       Use conda environments
        apptainer   Use Apptainer containers
        singularity Use Singularity containers
        docker      Use Docker containers
        test        Run with bundled test data

    Use --help to show this message.
    """.stripIndent()
    System.exit(0)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { METASSEMBLE } from './subworkflows/metassembly/main'
include { POLISH } from './subworkflows/polish/main'
include { EMAIL_RESULTS } from './modules/custom/email/main'
include { TWOBIT_TO_FA } from './modules/custom/ucsc/twobittofa/main'
include { CHROMSIZE } from './modules/custom/chromsize/main'
include { BED2GTF } from './modules/custom/bed2gtf/main'
include { GXF2BED } from './modules/custom/gxf2bed/main'
include { GUNZIP as GUNZIP_FASTA } from './modules/custom/gunzip/main'
include { GUNZIP as GUNZIP_GTF } from './modules/custom/gunzip/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def validateRun() {
    def errors = []
    if (!params.input_dir)   errors << "  --input_dir is required"
    if (!params.genome)      errors << "  --genome is required"
    if (!params.annotation)  errors << "  --annotation is required"

    if (errors) {
        log.error "Parameter validation failed:\n${errors.join('\n')}"
        System.exit(1)
    }
}

def validateFromPolishing() {
    def errors = []
    if (!params.polish_path)   errors << "  --polish_path is required"
    if (!params.genome)      errors << "  --genome is required"
    if (!params.annotation)  errors << "  --annotation is required"

    if (errors) {
        log.error "Parameter validation failed:\n${errors.join('\n')}"
        System.exit(1)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


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
        log.error "ERROR: Pipeline failed. Please refer to github issues: https://github.com/alejandrogzi/metassembly/issues"
    }

    workflow.onComplete {
        log.info "\nPipeline completed successfully!"
    }
}

// ── Default: full pipeline ─────────────────────────────────────────────────
workflow FULL_RUN {
    // validateRun()

    METASSEMBLE (
        params.input_dir,
        params.genome,
        params.annotation,
        params.star_index_path,
        params.star_ignore_gtf_for_index,
        params.deacon_index_path,
        params.deacon_download_index,
        params.deacon_make_single_index,
        params.star_make_coverage,
        params.skip_assembly,
        params.splice_scores_dir,
        params.output_dir
    )

    if (!params.skip_assembly) {
        POLISH (
            METASSEMBLE.out.metassembly,
            METASSEMBLE.out.genome,
            METASSEMBLE.out.annotation,
            params.repeats,
            params.splice_scores_dir,
        )
    }

    PIPELINE_COMPLETION (
        params.email_to,
        params.email_on_fail,
        params.plaintext_email,
        params.output_dir,
        params.use_mailx,
        METASSEMBLE.out.samplesheet
    )
}

// ── Checkpoint: start from polishing step (skip metassembly) ─────────────────────
workflow FROM_POLISHING {
    // validateFromPolishing()

    ch_versions = Channel.empty()

    def genome_file = file(params.genome, checkIfExists: true)
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
    if (params.annotation.endsWith('.gz') || parmas.annotation.endsWith('.gtf')) {
      Channel.value(file(params.annotation, checkIfExists: true))
        .map { it -> [ [ id: it.baseName ], it ] }
        .set { ch_gtf }

      if (params.annotation.endsWith('.gz')) {
         GUNZIP_GTF(ch_gtf)

         ch_gtf = GUNZIP_GTF.out.gunzip
         ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
      }

      // INFO: converting to bed
      GXF2BED(ch_gtf)
      ch_bed = GXF2BED.out.bed

      ch_versions = ch_versions.mix(GXF2BED.out.versions)
    } else if (params.annotation.endsWith('.bed')) {
      // INFO: converting to gtf
      Channel.value(file(params.annotation, checkIfExists: true))
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

    POLISH (
        params.polish_path,
        ch_fasta,
        ch_bed,
        params.repeats,
        params.splice_scores_dir,
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow XASM {
    if (params.from == "polish") {
        // ── Checkpoint: start from polishing step (skip metassembly) ─────────────────────
        log.info "Resuming from ${params.from} checkpoint — skipping LASTZ + chain building"
        FROM_POLISHING()
    } else {
        // ── Default: full pipeline ─────────────────────────────────────────────────
        log.info "Starting full pipeline — skipping checkpoints"
        FULL_RUN()
    }
}

workflow {
    XASM()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION HANDLER
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (workflow.success) {
        def results_dir = new File(params.output_dir as String, '10_results')
        def final_beds = results_dir.exists() ? (results_dir.listFiles()?.findAll { it.name.endsWith('.bed') } ?: []) : []
        log.info "Pipeline completed successfully!"
        if (final_beds) {
            log.info "Final predictions: ${final_beds.collect { it.toString() }.join(', ')}"
        } else {
            log.warn "Pipeline reported success but final bed file was not produced - check that all steps ran"
        }
        log.info "Run time   : ${workflow.duration}"
    } else {
        log.error "Pipeline FAILED — ${workflow.errorMessage}"
    }
}

workflow.onError {
    log.error "Pipeline error: ${workflow.errorMessage}"
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
