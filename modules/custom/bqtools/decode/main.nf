/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BQTOOLS_DECODE — Decode a columnar BINSEQ (CBQ) file back to gzipped FASTQ.
    A paired .cbq expands into ${prefix}_1.fastq.gz + ${prefix}_2.fastq.gz so STAR's
    even/odd mate split and --readFilesCommand zcat keep working. Single-end writes
    only ${prefix}_1.fastq.gz.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process BQTOOLS_DECODE {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' :
        'ghcr.io/hillerlab/bqtools:latest' }"

    input:
    tuple val(meta), path(cbq)
    val delete_input

    output:
    tuple val(meta), path("*_{1,2}.fastq.gz"), emit: fastq
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        bqtools \\
            decode \\
            -o ${prefix}_1.fastq.gz \\
            -f q \\
            -T $task.cpus \\
            $args \\
            $cbq

        if [ ${delete_input} == "true" ]; then
            for file in $cbq; do
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
            bqtools: \$(bqtools --version | sed 's/bqtools //g')
        END_VERSIONS
        """
    } else {
        """
        bqtools \\
            decode \\
            --prefix decode_tmp \\
            -f q \\
            -c g \\
            -T $task.cpus \\
            $args \\
            $cbq

        mv decode_tmp_R1.fq.gz ${prefix}_1.fastq.gz
        mv decode_tmp_R2.fq.gz ${prefix}_2.fastq.gz

        if [ ${delete_input} == "true" ]; then
            for file in $cbq; do
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
            bqtools: \$(bqtools --version | sed 's/bqtools //g')
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def stub_reads = meta.single_end
        ? "touch ${prefix}_1.fastq.gz"
        : "touch ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz"
    """
    ${stub_reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bqtools: \$(bqtools --version | sed 's/bqtools //g')
    END_VERSIONS
    """
}
