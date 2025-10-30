process FILTER_JUNCTIONS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(junctions)
    val discard_non_canonical_junctions
    val min_read_support
    val discard_canonical_junctions

    output:
    tuple val(meta), path("SJ_out_filtered.tab"), emit: filtered_junctions
    tuple val(meta), path("SJ_out_filtered.tab"), env(LINE_COUNT), emit: filtered_junctions_count
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Handle discard_non_canonical_junctions parameter
    def canonical_filter = discard_non_canonical_junctions ? '\$5 > 0' : '\$5 >= 0'

    // Handle min_read_support parameter
    def min_support = min_read_support ?: 0
    def read_support_filter = "\$7 > ${min_support}"

    // Handle discard_canonical_junctions parameter - only add if provided
    def motif_filter = discard_canonical_junctions != null && discard_canonical_junctions != '' ? " && \$6==${discard_canonical_junctions}" : ""

    // Construct the awk filter
    def awk_filter = "(${canonical_filter} && ${read_support_filter}${motif_filter})"

    """
    cat ${junctions} | awk '${awk_filter}' | cut -f1-6 | sort | uniq > SJ_out_filtered.tab

    LINE_COUNT=\$(wc -l < SJ_out_filtered.tab)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk --version | head -n1 | sed 's/GNU Awk //; s/,.*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch SJ_out_filtered.tab
    LINE_COUNT=0

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk --version | head -n1 | sed 's/GNU Awk //; s/,.*//')
    END_VERSIONS
    """
}
