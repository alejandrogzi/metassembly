/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FXSPLIT — Split FASTX/FASTQ reads into chunks for parallel processing.
    Partitions read files into smaller chunks to enable parallel processing
    across multiple CPUs. Set meta.headers or task.ext.headers to use
    fxsplit -H (one output file per FASTA record).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FXSPLIT {
    tag "$meta.id"
    label 'process_low'

    // ponytail: no conda env — the ANNEVO branch is container-only. Add an
    // environment.yml when someone needs -profile conda here.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' :
        'ghcr.io/alejandrogzi/fxsplit:latest' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("chunks/gz/*.gz")      , optional: true, emit: fastx_gz
    tuple val(meta), path("chunks/fx/*.f*")      , optional: true, emit: fastx
    path  "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args   ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def chunks        = task.ext.chunks ?: 400000
    def header_mode   = (meta.headers == true) || (task.ext.headers ? true : false)
    def mode_flag     = header_mode ? '-H' : "-c ${chunks}"
    def suffix_flag   = header_mode ? '' : "--suffix ${prefix}"
    def gzip          = reads.name.endsWith('.gz') ? true : false
    """
    fxsplit \\
        $args \\
        -f $reads \\
        ${mode_flag} \\
        -t $task.cpus \\
        -C \\
        ${suffix_flag}

    if [ ${meta.singleton} == true ]; then
        for f in chunks/*fasta.gz; do
            if [ -f "\$f" ]; then
                mv "\$f" "\${f%}.singleton.fasta.gz"
            fi
        done

        for f in chunks/*fastq.gz; do
            if [ -f "\$f" ]; then
                mv "\$f" "\${f%}.singleton.fastq.gz"
            fi
        done

        for f in chunks/*fasta; do
            if [ -f "\$f" ]; then
                mv "\$f" "\${f%}.singleton.fasta"
            fi
        done

        for f in chunks/*fastq; do
            if [ -f "\$f" ]; then
                mv "\$f" "\${f%}.singleton.fastq"
            fi
        done
    fi

    if [ $gzip == true ]; then
        mkdir chunks/gz
        mv chunks/*fast*.gz chunks/gz 2>/dev/null || mv chunks/*.gz chunks/gz
    else
        mkdir chunks/fx
        if [ "${header_mode}" = "true" ]; then
            find chunks -maxdepth 1 -type f -exec mv {} chunks/fx/ \\;
        else
            mv chunks/*fast* chunks/fx
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fxsplit: \$( fxsplit --version | head -n 1 | sed 's/fxsplit //g' | sed 's/ (.*//g' )
    END_VERSIONS
    """

    stub:
    // ponytail: the upstream stub touched `chunks/*fasta`, which creates a file
    // literally named "*fasta" that matches neither `chunks/fx/*.f*` nor
    // `chunks/gz/*.gz` — every downstream scatter then saw zero chunks. One
    // deterministic record (`chrStub`, the same name ANNEVO_MANIFEST's stub
    // reports) is enough to exercise the graph under -stub-run.
    def gzip   = reads.name.endsWith('.gz') ? true : false
    """
    mkdir chunks

    if [ $gzip == true ]; then
        mkdir chunks/gz
        touch chunks/gz/chrStub.fa.gz
    else
        mkdir chunks/fx
        touch chunks/fx/chrStub.fa
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fxsplit: \$( fxsplit --version | head -n 1 | sed 's/fxsplit //g' | sed 's/ (.*//g' )
    END_VERSIONS
    """
}
