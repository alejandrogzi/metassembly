/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { STAR_GENOMEGENERATE } from '../../modules/nf-core/star/genomegenerate/main'
include { UNTAR as UNTAR_STAR_INDEX } from '../../modules/nf-core/untar/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREPARE_GENOME_STAR {
    take:
        fasta                // file: /path/to/genome.fasta
        gtf                  // file: /path/to/genome.gtf
        star_index_path      // path: /path/to/star/index/
        star_ignore_sjdbgtf  // val: boolean

    main:
        // Versions collector + init
        ch_versions = Channel.empty()
        ch_gtf = Channel.of([])

        ch_fasta = Channel.value(file(fasta, checkIfExists: true))

        if (!star_ignore_sjdbgtf) {
            ch_gtf = Channel.value(file(gtf, checkIfExists: true))
        }

        ch_star_index = Channel.empty()
        if (star_index_path) {
            if (star_index_path.endsWith('.tar.gz')) {
                ch_star_index = UNTAR_STAR_INDEX([[:], file(star_index_path, checkIfExists: true)]).untar.map { it[1] }
                ch_versions = ch_versions.mix(UNTAR_STAR_INDEX.out.versions)
            } else {
                ch_star_index = Channel.value(file(star_index_path, checkIfExists: true))
            }
        } else {
            ch_star_index = STAR_GENOMEGENERATE(
                ch_fasta.map { [[:], it] },
                ch_gtf.map { [[:], it] }
            ).index.map { it[1] }

            ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
        }

    emit:
        star_index = ch_star_index // channel: path(star/index)
        star_gtf = ch_gtf // channel: path(star/sjdb.gtf)
        versions = ch_versions // channel: [ versions.yml ]
}
