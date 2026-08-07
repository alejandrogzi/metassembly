// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

process DEACON_FILTER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' :
        'ghcr.io/hillerlab/deacon-cbq@sha256:4d684c770b08551fdc89588598c87f3f038712388b2e56eabbc94b18f2a35844' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta1), path(index)

    output:
    tuple val(meta), path("*.deacon.{fastq.gz,cbq}"), emit: reads
    tuple val(meta), path("*.deacon.log")           , emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // A .cbq holds both mates in one file and deacon rejects --output2 for CBQ output;
    // the .cbq suffix on --output is what selects CBQ encoding.
    def is_cbq = (reads instanceof List ? reads[0] : reads).name.endsWith('.cbq')
    def outputs = is_cbq
        ? "--output ${prefix}.deacon.cbq"
        : (meta.single_end
            ? "--output ${prefix}_1.deacon.fastq.gz"
            : "--output ${prefix}_1.deacon.fastq.gz --output2 ${prefix}_2.deacon.fastq.gz")
    """
    deacon \\
        filter \\
        --threads $task.cpus \\
        $outputs \\
        $args \\
        $index \\
        $reads \\
        > ${prefix}.deacon.log 2>&1

    if [ ${params.deacon_keep_fastp_fastq} == false ]; then
        # Resolve symlinks and delete actual files
        for file in $reads; do
            if [ -L "\$file" ]; then
                realpath=\$(readlink -f "\$file")
                rm -f "\$realpath"
            else
                rm -f "\$file"
            fi
        done
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deacon: \$(deacon --version | head -n1 | sed 's/deacon //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def stub_reads = (reads instanceof List ? reads[0] : reads).name.endsWith('.cbq')
        ? "touch ${prefix}.deacon.cbq"
        : "touch ${prefix}_1.deacon.fastq.gz ${prefix}_2.deacon.fastq.gz"
    """
    touch ${prefix}.idx
    ${stub_reads}
    touch ${prefix}.deacon.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deacon: \$(deacon --version | head -n1 | sed 's/deacon //g')
    END_VERSIONS
    """
}
