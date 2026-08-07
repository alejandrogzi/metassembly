/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FQ_LINT as FQ_LINT_AT_START } from '../../modules/nf-core/fq/lint/main'
include { FQ_LINT as FQ_LINT_AFTER_TRIMMING } from '../../modules/nf-core/fq/lint/main'
include { FQ_LINT as FQ_LINT_AFTER_DECONTAMINATION } from '../../modules/nf-core/fq/lint/main'
include { FASTP } from '../../modules/nf-core/fastp/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BQC } from '../../modules/custom/bqc/workflow/main'
include { DEACON_FILTER } from '../../modules/custom/deacon/filter/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREPROCESS_READS {
    take:
        fastqs
        deacon_index

    main:
        // Versions collector + init
        ch_versions = Channel.empty()
        ch_linting_logs = Channel.empty()

        // fastp cannot read CBQ and bqc only reads CBQ, so the QC tool is derived
        // from the input format rather than configured.
        fastqs
            .branch { _meta, reads ->
                cbq: isCbq(reads)
                fastq: true
            }
            .set { ch_by_format }

        if (!params.fq_skip_linting_at_start) {
            FQ_LINT_AT_START(
                ch_by_format.fastq
            )

            ch_versions = ch_versions.mix(FQ_LINT_AT_START.out.versions)
            ch_linting_logs = ch_linting_logs.mix(FQ_LINT_AT_START.out.lint)
        }

        FASTP(
            ch_by_format.fastq,
            [], // adapter_fasta
            params.fastp_discard_trimmed_pass,
            params.fastp_save_trimmed_fail,
            params.fastp_save_merged
        )

        trim_json = FASTP.out.json
        trim_html = FASTP.out.html
        trim_log = FASTP.out.log

        ch_versions = ch_versions.mix(FASTP.out.versions.first())

        BQC(
            ch_by_format.cbq
        )

        ch_versions = ch_versions.mix(BQC.out.versions.first())

        def minTrimmedReads = params.fastp_min_trimmed_reads.toLong()

        FASTP.out.reads.join(trim_json)
            .map { meta, _reads, json -> [meta, _reads, getFastpReadsAfterFiltering(json, minTrimmedReads)] }
            .mix(
                BQC.out.reads.join(BQC.out.report)
                    .map { meta, _reads, json -> [meta, _reads, getBqcReadsAfterFiltering(json, minTrimmedReads)] }
            )
            .set { ch_num_trimmed_reads }

        trim_json
            .map { meta, json -> [meta, getFastpReadsAfterFilteringAsPercentage(json, minTrimmedReads)] }
            .mix(
                BQC.out.report
                    .map { meta, json -> [meta, getBqcReadsAfterFilteringAsPercentage(json, minTrimmedReads)] }
            )
            .set { ch_num_trimmed_reads_percent }

        ch_num_trimmed_reads
            .filter { meta, _reads, num_reads ->
                def keepSample = num_reads >= minTrimmedReads
                if (!keepSample) {
                    def sampleId = meta?.id ?: meta
                    log.warn "[PREPROCESS_READS] Discarding sample ${sampleId} after trimming: ${num_reads} reads < min_trimmed_reads (${minTrimmedReads})"
                }
                return keepSample
            }
            .map { meta, _reads, _num_reads -> [meta, _reads] }
            .set { ch_trimmed_reads }

        ch_num_trimmed_reads
            .map { meta, _reads, num_reads -> [meta, num_reads] }
            .set { trim_read_count }

        trim_json
            .map { meta, json -> [meta, getFastpAdapterSequence(json)] }
            .set { ch_adapter_seq }

        if (!params.fq_skip_linting_after_trimming) {
            // fq lint is FASTQ-only; CBQ integrity is covered by bqc/bqtools
            FQ_LINT_AFTER_TRIMMING(
                ch_trimmed_reads.filter { _meta, reads -> !isCbq(reads) }
            )

            ch_versions = ch_versions.mix(FQ_LINT_AFTER_TRIMMING.out.versions)
            ch_linting_logs = ch_linting_logs.mix(FQ_LINT_AFTER_TRIMMING.out.lint)
        }

        ch_deacon_out = DEACON_FILTER(
            ch_trimmed_reads,
            deacon_index
        )

        ch_versions = ch_versions.mix(DEACON_FILTER.out.versions)

        ch_deacon_log = DEACON_FILTER.out.log
        ch_deacon_log
            .map { meta, log -> [meta, getDeaconRetainedPercent(log)] }
            .set { ch_deacon_discarded_seqs }

        if (!params.fq_skip_linting_after_deacon) {
            FQ_LINT_AFTER_DECONTAMINATION(
                ch_deacon_out.reads.filter { _meta, reads -> !isCbq(reads) }
            )

            ch_versions = ch_versions.mix(FQ_LINT_AFTER_DECONTAMINATION.out.versions)
            ch_linting_logs = ch_linting_logs.mix(FQ_LINT_AFTER_DECONTAMINATION.out.lint)
        }

    emit:
        processed_reads = ch_deacon_out.reads
        adapter_seq = ch_adapter_seq
        lint_logs = ch_linting_logs
        num_trimmed_reads = trim_read_count
        num_trimmed_reads_percent = ch_num_trimmed_reads_percent
        deacon_discarded_seqs = ch_deacon_discarded_seqs
        fastp_html = FASTP.out.html
        fastp_json = FASTP.out.json
        versions = ch_versions
}

//
// Function that parses and returns the number of reads after filtering from FastP log output
//
def getFastpReadsAfterFiltering(json_file, min_num_reads) {
    if (workflow.stubRun) {
        return min_num_reads
    }

    def json = new groovy.json.JsonSlurper().parseText(json_file.text).get('summary') as Map
    return json['after_filtering']['total_reads'].toLong()
}


//
// Function that parses and returns the number of reads after filtering from FastP log output
//
def getFastpReadsAfterFilteringAsPercentage(json_file, min_num_reads) {
    if (workflow.stubRun) {
        return min_num_reads
    }

    def json = new groovy.json.JsonSlurper().parseText(json_file.text).get('summary') as Map

    def total_reads_before_filtering = json['before_filtering']['total_reads'].toLong()
    def total_reads_after_filtering = json['after_filtering']['total_reads'].toLong()

    return (total_reads_after_filtering / total_reads_before_filtering) * 100
}

//
// Function that reports whether a reads entry is CBQ rather than FASTQ.
// Mirrors metassembly/main.nf: the format comes from the file, not from meta.
//
def isCbq(reads) {
    def first = reads instanceof List ? reads[0] : reads
    return first.name.endsWith('.cbq')
}

//
// Function that parses and returns the number of reads after filtering from the bqc report
//
// NOTE: bqc counts *records* (a paired record carries both mates) while fastp counts
// *reads*, so a paired file is doubled to keep both branches on the same scale --
// otherwise fastp_min_trimmed_reads would discard CBQ samples at twice the threshold.
//
def getBqcReadsAfterFiltering(json_file, min_num_reads) {
    if (workflow.stubRun) {
        return min_num_reads
    }

    def counts = new groovy.json.JsonSlurper().parseText(json_file.text).get('counts') as Map
    def mates = (counts['r2_bases_in'] as Long) > 0 ? 2 : 1

    return (counts['records_out'] as Long) * mates
}

//
// Function that parses and returns the percentage of reads after filtering from the bqc report
//
def getBqcReadsAfterFilteringAsPercentage(json_file, min_num_reads) {
    if (workflow.stubRun) {
        return min_num_reads
    }

    def counts = new groovy.json.JsonSlurper().parseText(json_file.text).get('counts') as Map

    def records_before_filtering = counts['records_in'] as Long
    def records_after_filtering = counts['records_out'] as Long

    // Ratio is mate-count independent, so no doubling is needed here
    return records_before_filtering ? (records_after_filtering / records_before_filtering) * 100 : 0
}

//
// Function that parses and returns the adapter sequence from FastP log output
//
def getFastpAdapterSequence(json_file) {
    // Handle stub runs
    if (workflow.stubRun) {
        return ""
    }

    def json = new groovy.json.JsonSlurper().parseText(json_file.text) as Map
    try {
        return json['adapter_cutting']['read1_adapter_sequence']
    } catch (Exception ex) {
        return ""
    }
}

//
// Function that parses and returns the retained sequence percentage from Deacon log output
//
def getDeaconRetainedPercent(deacon_log) {
    def retained_percent = 0
    def pattern = /Retained\s+\d+\/\d+\s+sequences\s+\(([\d\.]+)%\)/

    deacon_log.eachLine { line ->
        def matcher = line =~ pattern
        if (matcher) {
            retained_percent = matcher[0][1].toFloat()
        }
    }

    return retained_percent
}
