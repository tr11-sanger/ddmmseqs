process DIAMOND_MAKEDB_FROM_PIPE {
    tag "${meta.id}"
    label 'process_medium'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/diamond:2.1.12--hdb4b4cc_1'
        : 'biocontainers/diamond:2.1.12--hdb4b4cc_1'}"

    input:
    tuple val(meta), path(filelist)

    output:
    tuple val(meta), path("*.dmnd"), emit: db
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    cat ${filelist} \\
    | while read fp; do
        if [[ \$fp =~ \\.gz\$ ]]; then
            gunzip -c \$fp
        else
            cat \$fp
        fi
    done \\
    | diamond makedb \\
        --threads ${task.cpus} \\
        -d ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | tail -n 1 | sed 's/^diamond version //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "${args}"
    touch ${prefix}.dmnd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | tail -n 1 | sed 's/^diamond version //')
    END_VERSIONS
    """
}
