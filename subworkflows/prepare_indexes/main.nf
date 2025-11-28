/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPARE_GENOME_STAR } from '../star_index/main'
include { PREPARE_DEACON_INDEX } from '../deacon_index/main'
include { TWOBIT_TO_FA } from '../../modules/custom/twobit/main'
include { GUNZIP_FASTA } from '../../modules/custom/gunzip/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREPARE_INDEXES {
    take:
        fasta                               // file: /path/to/genome.fasta
        // star_index
        gtf                                 // file: /path/to/genome.gtf
        star_index_path                     // path: /path/to/star/index/
        star_ignore_gtf_for_index           // val: boolean
        // deacon_index
        index_path                          // val: path(deacon/index)
        download_index                      // val: boolean
        make_single_index                   // val: boolean
        multi_index_additional_genome_paths // val: list(path(genome))

    main:
        ch_versions = Channel.empty()
        ch_fasta_preprocessed = Channel.empty()

        def fasta_file = file(fasta, checkIfExists: true)
        def fasta_path = fasta_file.toString()

        // INFO: if fasta is .2bit or .gz, convert or uncompress it
        if (fasta_path.endsWith(".2bit")) {
            TWOBIT_TO_FA([[:], fasta_file])
            ch_fasta_preprocessed = TWOBIT_TO_FA.out.fasta.map { it[1] }
            ch_versions = ch_versions.mix(TWOBIT_TO_FA.out.versions)
        } else if (fasta_path.endsWith(".gz")) {
            GUNZIP_FASTA([[:], fasta_file])
            ch_fasta_preprocessed = GUNZIP_FASTA.out.fasta.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
        } else {
            ch_fasta_preprocessed = Channel.value(fasta_file)
        }

        PREPARE_GENOME_STAR(
            ch_fasta_preprocessed,
            gtf,
            star_index_path,
            star_ignore_gtf_for_index,
        )

        ch_versions = ch_versions.mix(PREPARE_GENOME_STAR.out.versions)

        PREPARE_DEACON_INDEX(
            index_path,
            download_index,
            make_single_index,
            multi_index_additional_genome_paths,
            ch_fasta_preprocessed
        )

        ch_versions = ch_versions.mix(PREPARE_DEACON_INDEX.out.versions)

    emit:
        star_index = PREPARE_GENOME_STAR.out.star_index
        star_gtf = PREPARE_GENOME_STAR.out.star_gtf
        deacon_index = PREPARE_DEACON_INDEX.out.deacon_index
        versions = ch_versions
}
