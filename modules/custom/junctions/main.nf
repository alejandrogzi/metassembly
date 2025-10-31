process JOIN_JUNCTIONS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10-slim' :
        'biocontainers/python:3.10-slim' }"

    input:
    tuple val(meta), path(junctions)
    val min_junction_len
    val min_junction_coverage

    output:
    tuple val(meta), path("*.tab"), emit: filtered_junctions
    tuple val(meta), path("*.tab"), env(LINE_COUNT), emit: filtered_junctions_count
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def junction_files = junctions.join(' ')

    """
    join_junctions.py -j ${junction_files} -l ${min_junction_len} -m ${min_junction_coverage} -o .
    LINE_COUNT=\$(wc -l < ALL_SJ_out_filtered.tab)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        join_junctions.py: \$(join_junctions.py --version | sed 's/join_junctions.py //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ALL_SJ_out_filtered.tab
    LINE_COUNT=0

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        join_junctions.py: \$(join_junctions.py --version | sed 's/join_junctions.py //')
    END_VERSIONS
    """
}
