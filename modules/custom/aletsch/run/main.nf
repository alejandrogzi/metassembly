process ALETSCH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/aletsch:1.1.3--h9f5acd7_0' :
        'biocontainers/aletsch:1.1.3--h9f5acd7_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.gtf")       , emit: gtf
    tuple val(meta), path("*profile")    , emit: profile
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def library_type = meta.strandedness ?: 'unstranded'

    """
    # Create necessary directories
    mkdir -p ${prefix}_profile
    mkdir -p ${prefix}_gtf

    # Create sample info file
    echo "${bam}\t${bai}\t${library_type}" > ${prefix}.info

    # Run Aletsch profile generation
    aletsch \\
        --profile \\
        -i ${prefix}.info \\
        -p ${prefix}_profile \\
        $args

    # Run Aletsch assembly
    aletsch \\
        -i ${prefix}.info \\
        -o ${prefix}_gtf/${prefix}.gtf \\
        -p ${prefix}_profile \\
        -d ${prefix}_gtf \\
        $args

    # Move output to current directory
    mv ${prefix}_gtf/${prefix}.gtf ${prefix}.gtf
    mv ${prefix}_profile ${prefix}.profile

    rename_gtf.py \\
        -g ${prefix}.gtf \\
        -p ${prefix} \\
        -o ${prefix}.renamed.gtf

    rm ${prefix}.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aletsch: \$(aletsch --version 2>&1 | sed 's/^.*aletsch //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.renamed.gtf

    mkdir -p ${prefix}.profile
    touch ${prefix}.profile/0.profile

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aletsch: \$(aletsch --version 2>&1 | sed 's/^.*aletsch //; s/ .*\$//')
    END_VERSIONS
    """
}
