/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPARE_GENOME_STAR } from '../star_index/main'
include { PREPARE_DEACON_INDEX } from '../deacon_index/main'
include { TWOBIT_TO_FA } from '../../modules/custom/ucsc/twobittofa/main'
include { GUNZIP as GUNZIP_FASTA } from '../../modules/custom/gunzip/main'
include { CHROMSIZE } from '../../modules/custom/chromsize/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREPARE_INDEXES {
    take:
        genome                               // file: /path/to/genome.{2bit/fasta}
        // star_index
        gtf                                 // file: /path/to/genome.gtf
        star_index_path                     // path: /path/to/star/index/
        star_ignore_gtf_for_index           // val: boolean
        // deacon_index
        index_path                          // val: path(deacon/index)
        download_index                      // val: boolean
        make_single_index                   // val: boolean

    main:
        ch_versions = Channel.empty()

        def genome_file = file(genome, checkIfExists: true)
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

        PREPARE_GENOME_STAR(
            ch_fasta,
            gtf,
            star_index_path,
            star_ignore_gtf_for_index,
        )

        ch_versions = ch_versions.mix(PREPARE_GENOME_STAR.out.versions)

        PREPARE_DEACON_INDEX(
            index_path,
            download_index,
            make_single_index,
            ch_fasta
        )

        ch_versions = ch_versions.mix(PREPARE_DEACON_INDEX.out.versions)

    emit:
        star_index = PREPARE_GENOME_STAR.out.star_index
        star_gtf = PREPARE_GENOME_STAR.out.star_gtf
        deacon_index = PREPARE_DEACON_INDEX.out.deacon_index
        chrom_sizes = ch_chrom_sizes
        versions = ch_versions
}
