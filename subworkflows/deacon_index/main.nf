/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DEACON_INDEX } from '../../modules/nf-core/deacon/index/main'
include { DEACON_MULTI_INDEX } from '../../modules/custom/deacon/multindex/main'
include { DEACON_MULTI_INDEX as DEACON_INDEX_WITH_BACKGROUND } from '../../modules/custom/deacon/multindex/main'
include { WGET } from '../../modules/nf-core/wget/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREPARE_DEACON_INDEX {
    take:
        index_path
        download_index
        make_single_index
        multi_index_additional_genome_paths
        fasta
    main:
        ch_versions = Channel.empty()
        ch_deacon_index = Channel.empty()

        ch_fasta = (fasta instanceof groovyx.gpars.dataflow.DataflowReadChannel) ?
            fasta.map { file(it, checkIfExists: true) } : Channel.value(file(fasta, checkIfExists: true))

        def additional_genomes = multi_index_additional_genome_paths ?: []

        if (index_path) {
            if (!download_index) {
                // INFO: return index_path -> user has provided index
                ch_deacon_index = Channel.value(file(index_path, checkIfExists: true))
            } else {
                // INFO: download index from index_path, extract file from tuple
                ch_deacon_index = WGET(
                    [
                        meta: ["id": file(index_path).name],
                        path: index_path
                    ]
                ).outfile.map { meta, file -> file }  // Extract file from [meta, file] tuple
                ch_versions = ch_versions.mix(WGET.out.versions)
            }
        } else {
            if (make_single_index) {
                // INFO: create index using fasta
                if (!params.deacon_single_index_use_background) {
                    def deacon_output = DEACON_INDEX(
                        ch_fasta
                    )
                    ch_deacon_index = deacon_output.index
                    ch_versions = ch_versions.mix(deacon_output.versions)
                } else {
                    def deacon_output = DEACON_INDEX_WITH_BACKGROUND(
                        ch_fasta
                    )
                    ch_deacon_index = deacon_output.index
                    ch_versions = ch_versions.mix(deacon_output.versions)
                }
            } else {
                ch_fasta
                    .map { fasta_file ->
                        def genomes = [fasta_file]
                        genomes.addAll(additional_genomes.collect { file(it, checkIfExists: true) })
                        genomes.join(',')
                    }
                    .set { ch_multi_index_genome_paths }
                // INFO: create multi-index using multi_index_genome_paths + fasta
                def deacon_output = DEACON_MULTI_INDEX(
                    ch_multi_index_genome_paths
                )
                ch_deacon_index = deacon_output.index.map { meta, file -> file }
                ch_versions = ch_versions.mix(deacon_output.versions)
            }
        }
    emit:
        deacon_index = ch_deacon_index // channel: path(deacon/index)
        versions = ch_versions // channel: [ versions.yml ]
}
