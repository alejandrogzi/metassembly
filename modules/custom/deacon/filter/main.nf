process DEACON_FILTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deacon:0.12.0--h4349ce8_0':
        'biocontainers/deacon:0.12.0--h4349ce8_0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta1), path(index)

    output:
    tuple val(meta), path("*.deacon.fastq.gz"), emit: reads
    tuple val(meta), path("*.deacon.log")     , emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def out_fq1 = "--output ${prefix}_1.deacon.fastq.gz"
    def out_fq2 = "--output2 ${prefix}_2.deacon.fastq.gz"
    """
    deacon \\
        filter \\
        --threads $task.cpus \\
        $out_fq1 \\
        $out_fq2 \\
        $args \\
        $index \\
        $reads \\
        > ${prefix}.deacon.log 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deacon: \$(deacon --version | head -n1 | sed 's/deacon //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.idx
    touch ${prefix}_1.deacon.fastq.gz
    touch ${prefix}_2.deacon.fastq.gz
    touch ${prefix}.deacon.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deacon: \$(deacon --version | head -n1 | sed 's/deacon //g')
    END_VERSIONS
    """
}
