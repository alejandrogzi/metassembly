/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TRANSMETA — Multi-sample RNA-seq transcript meta-assembly.
    Simultaneously assembles RNA-seq reads of multiple samples into a unified set
    of transcripts (GTF) and a set of transcripts for each individual sample.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process TRANSMETA {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' :
        'ghcr.io/hillerlab/transmeta@sha256:d0a7b62ec69620f6ddc31c1b753f0b08c046241473c4eabb13d29f6384d00907' }"

    input:
    tuple val(meta), path(bams, stageAs: "bams/*")
    tuple val(meta1), path(annotation)

    output:
    tuple val(meta), path("transmeta_outdir/TransMeta.gtf"),            emit: gtf
    tuple val(meta), path("transmeta_outdir/TransMeta-[0-9]*.gtf"),     emit: meta_gtf
    tuple val(meta), path("transmeta_outdir/sample_gtfs/*.gtf"),        emit: sample_gtf
    tuple val(meta), path("transmeta_outdir/TransMeta-AG*.gtf"),        optional: true, emit: ag_gtf
    tuple val(meta), path("transmeta_outdir/counts.tsv"),               emit: counts
    path "versions.yml",                                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def strand = meta.single_end ?
        (meta.strandedness == 'reverse' ? 'single_reverse' : meta.strandedness == 'forward' ? 'single_forward' : 'single_unstranded') :
        (meta.strandedness == 'reverse' ? 'second' : meta.strandedness == 'forward' ? 'first' : 'unstranded')

    def annotation_arg = params.transmeta_use_annotation && annotation ? "-g $annotation" : ''

    """
    ls bams/*.bam > bam.list

    TransMeta \\
        -B bam.list \\
        -s $strand \\
        -o transmeta_outdir \\
        -p $task.cpus \\
        $annotation_arg \\
        $args

    # Map the per-sample assemblies back to sample names and count transcripts
    mkdir -p transmeta_outdir/sample_gtfs
    > transmeta_outdir/counts.tsv
    idx=1
    for gtf in transmeta_outdir/TransMeta.bam*.gtf; do
        [ -f "\$gtf" ] || continue
        name=\$(basename "\$(sed -n "\${idx}p" bam.list)")
        name=\${name%.bam}
        mv "\$gtf" "transmeta_outdir/sample_gtfs/\${name}.gtf"
        count=\$(grep -w 'transcript' "transmeta_outdir/sample_gtfs/\${name}.gtf" | wc -l)
        printf "%s\\t%s\\n" "\$name" "\$count" >> transmeta_outdir/counts.tsv
        idx=\$((idx+1))
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transmeta: \$(TransMeta -v 2>&1 | sed 's/^.*v[.]//; s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p transmeta_outdir/sample_gtfs
    touch transmeta_outdir/TransMeta.gtf
    touch transmeta_outdir/TransMeta-*.gtf
    touch transmeta_outdir/TransMeta-AG*.gtf
    touch transmeta_outdir/counts.tsv
    touch transmeta_outdir/sample_gtfs/*.gtf

    > transmeta_outdir/counts.tsv
    for bam in bams/*.bam; do
        name=\$(basename "\$bam")
        name=\${name%.bam}
        touch "transmeta_outdir/sample_gtfs/\${name}.gtf"
        printf "%s\\t%s\\n" "\$name" "0" >> transmeta_outdir/counts.tsv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transmeta: \$(TransMeta -v 2>&1 | sed 's/^.*v[.]//; s/ .*//')
    END_VERSIONS
    """
}
