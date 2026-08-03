// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

process RENAME_GTF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10.2' :
        'biocontainers/python:3.10.2' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.gtf"), emit: gtf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    rename_gtf.py \\
        -g ${gtf} \\
        -p ${prefix} \\
        -o ${prefix}.renamed.gtf

    cat <<-END_VERSIONS > versions.yml  
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        rename_gtf.py: \$(rename_gtf.py --version | sed 's/rename_gtf.py //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.renamed.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        rename_gtf.py: \$(rename_gtf.py --version | sed 's/rename_gtf.py //')
    END_VERSIONS
    """
}
