include { PREPARE_GENOME_STAR } from '../star_index/main'
include { PREPARE_DEACON_INDEX } from '../deacon_index/main'

workflow PREPARE_INDEXES {
    take:
        fasta // file: /path/to/genome.fasta
        // star_index
        gtf // file: /path/to/genome.gtf
        star_ignore_sjdbgtf // val: boolean
        star_index_path // path: /path/to/star/index/
        // deacon_index
        index_path // val: path(deacon/index)
        download_index // val: boolean
        make_single_index // val: boolean
        multi_index_additional_genome_paths // val: list(path(genome))

    main:
        ch_versions = Channel.empty()

        PREPARE_GENOME_STAR(
            fasta,
            gtf,
            star_ignore_sjdbgtf,
            star_index_path
        )

        ch_versions = ch_versions.mix(PREPARE_GENOME_STAR.out.versions)

        PREPARE_DEACON_INDEX(
            index_path,
            download_index,
            make_single_index,
            multi_index_additional_genome_paths,
            fasta
        )

        ch_versions = ch_versions.mix(PREPARE_DEACON_INDEX.out.versions)

    emit:
        star_index = PREPARE_GENOME_STAR.out.star_index
        star_gtf = PREPARE_GENOME_STAR.out.star_gtf
        deacon_index = PREPARE_DEACON_INDEX.out.deacon_index
        versions = ch_versions
}
