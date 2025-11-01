process BEAVER {
    tag "meta_assembly"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    path gtfs

    output:
    path "beaver_output/*gtf"        , emit: gtf
    path "beaver_output/*csv"        , emit: csv
    path "gtf_list.txt"              , emit: gtf_list
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "beaver_output"

    """
    # Create GTF list file
    for gtf in ${gtfs}; do
        echo "\${gtf}" >> gtf_list.txt
    done

    # Create output directory
    mkdir -p beaver_output

    # Run Beaver
    beaver \\
        gtf_list.txt \\
        ${prefix} \\
        $args \\
        > beaver.log 2>&1

    mv ${prefix}.gtf beaver_output/
    mv ${prefix}_feature.csv beaver_output/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beaver: \$(beaver --version 2>&1 | sed 's/^.*beaver //; s/ .*\$//' || echo "0.0.1")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "beaver_output"
    """
    mkdir -p beaver_output
    touch gtf_list.txt
    touch beaver_output/${prefix}.gtf
    touch beaver_output/${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beaver: 0.0.1
    END_VERSIONS
    """
}
