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

        // FQ_LINT_AT_START(
        //     fastqs
        // )

        // ch_versions = ch_versions.mix(FQ_LINT_AT_START.out.versions)
        // ch_linting_logs = ch_linting_logs.mix(FQ_LINT_AT_START.out.lint)

        FASTP(
            fastqs,
            [], // adapter_fasta
            params.fastp_discard_trimmed_pass,
            params.fastp_save_trimmed_fail,
            params.fastp_save_merged
        )

        trim_json = FASTP.out.json
        trim_html = FASTP.out.html
        trim_log = FASTP.out.log

        ch_versions = ch_versions.mix(FASTP.out.versions.first())

        FASTP.out.reads.join(trim_json)
            .map { meta, _reads, json -> [meta, _reads, getFastpReadsAfterFiltering(json, params.fastp_min_trimmed_reads.toLong())] }
            .set { ch_num_trimmed_reads }

        ch_num_trimmed_reads
            .filter { _meta, _reads, num_reads -> num_reads >= params.fastp_min_trimmed_reads.toLong() }
            .map { meta, _reads, _num_reads -> [meta, _reads] }
            .set { ch_trimmed_reads }

        ch_num_trimmed_reads
            .map { meta, _reads, num_reads -> [meta, num_reads] }
            .set { trim_read_count }

        trim_json
            .map { meta, json -> [meta, getFastpAdapterSequence(json)] }
            .set { ch_adapter_seq }

        // FQ_LINT_AFTER_TRIMMING(
        //     ch_trimmed_reads
        // )

        // ch_versions = ch_versions.mix(FQ_LINT_AFTER_TRIMMING.out.versions)
        // ch_linting_logs = ch_linting_logs.mix(FQ_LINT_AFTER_TRIMMING.out.lint)

        ch_deacon_out = DEACON_FILTER(
            ch_trimmed_reads,
            deacon_index.map { i -> [ [id: "deacon_index"], i ] }
        )

        ch_versions = ch_versions.mix(DEACON_FILTER.out.versions)

        ch_deacon_log = DEACON_FILTER.out.log
        ch_deacon_log
            .map { meta, log -> [meta, getDeaconRetainedPercent(log)] }
            .set { ch_deacon_discarded_seqs }

        // FQ_LINT_AFTER_DECONTAMINATION(
        //     ch_deacon_out.reads
        // )

        // ch_versions = ch_versions.mix(FQ_LINT_AFTER_DECONTAMINATION.out.versions)
        // ch_linting_logs = ch_linting_logs.mix(FQ_LINT_AFTER_DECONTAMINATION.out.lint)

    emit:
        processed_reads = ch_deacon_out.reads
        adapter_seq = ch_adapter_seq
        lint_logs = ch_linting_logs
        num_trimmed_reads = trim_read_count
        deacon_discarded_seqs = ch_deacon_discarded_seqs
        versions = ch_versions
}

def getFastpReadsAfterFiltering(json_file, min_num_reads) {
    if (workflow.stubRun) {
        return min_num_reads
    }

    def json = new groovy.json.JsonSlurper().parseText(json_file.text).get('summary') as Map
    return json['after_filtering']['total_reads'].toLong()
}

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
def getDeaconRetainedPercent(params, deacon_log) {
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
