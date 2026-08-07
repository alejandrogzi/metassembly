// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

process REMOVE_BAMS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    // NOTE: was 'biocontainers/bash:5.2.15', which does not exist on the configured
    // quay.io default registry. This process is gated on
    // !aletsch_keep_bam && !star_make_coverage, so the broken coordinate was never
    // exercised by the default (coverage-on) configuration.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'quay.io/nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(bams)

    output:
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    for file in ${bams}; do
        if [ -L "\$file" ]; then
            realpath=\$(readlink -f "\$file")
            rm -f "\$realpath"
        fi
        rm -f "\$file"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | sed 's/.* version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: 5.2.15
    END_VERSIONS
    """
}