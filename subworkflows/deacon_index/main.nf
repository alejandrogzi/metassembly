/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DEACON_INDEX } from '../../modules/nf-core/deacon/index/main'
include { DEACON_MULTI_INDEX } from '../../modules/custom/deacon/multindex/main'
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
        // Versions collector + init
        ch_versions = Channel.empty()
        ch_deacon_index = Channel.empty()

        if (index_path) {
            if (!download_index) {
                // INFO: return index_path -> user has provided index
                ch_deacon_index = Channel.value(file(index_path, checkIfExists: true))
            } else {
                // INFO: download index from index_path, emit multi-channel
                ch_deacon_index = WGET(
                    [
                        meta: ["id": file(index_path).name],
                        path: index_path
                    ]
                ).outfile
            }
        } else {
            if (make_single_index) {
                // INFO: create index using fasta
                ch_deacon_index = DEACON_INDEX(
                    fasta
                )
            } else {
                def genomes = []
                genomes << fasta
                genomes = genomes.concat(multi_index_additional_genome_paths)

                multi_index_genome_paths = genomes.collect { file(it, checkIfExists: true) }

                // INFO: make it 1 string separate by comma
                multi_index_genome_paths = multi_index_genome_paths.join(',').toString()

                // INFO: create multi-index using multi_index_genome_paths + fasta
                ch_deacon_index = DEACON_MULTI_INDEX(
                    multi_index_genome_paths
                )
            }
        }

    emit:
        deacon_index = ch_deacon_index // channel: path(deacon/index)
        versions = ch_versions // channel: [ versions.yml ]
}
