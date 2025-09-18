process MST_TO_CLUSTER {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    //     ? 'https://depot.galaxyproject.org/singularity/diamond:2.1.12--hdb4b4cc_1'
    //     : 'biocontainers/diamond:2.1.12--hdb4b4cc_1'}"

    input:
    tuple val(meta), path(filelists), path(node_linkage_lists)
    val target_cluster_size

    output:
    tuple val(meta), path('*_nodes.{tsv,tsv.gz}'), emit: nodes
    tuple val(meta), path('*_clusters.{tsv,tsv.gz}'), emit: clusters
    tuple val(meta), path('*_cluster_seqs/*.{faa,faa.gz}'), emit: cluster_seqs
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python ${moduleDir}/bin/mst_to_cluster.py \\
        --linkagelists ${node_linkage_lists} \\
        --filelists ${filelists} \\
        --out_nodes ${prefix}_nodes.tsv.gz \\
        --out_clusters ${prefix}_clusters.tsv.gz \\
        --out_cluster_seqs_dir ${prefix}_cluster_seqs \\
        --target_cluster_size ${target_cluster_size} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | tail -n 1 | sed 's/^diamond version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_nodes.txt.gz
    touch ${prefix}_clusters.tsv.gz
    mkdir ${prefix}_cluster_seqs
    touch ${prefix}_cluster_seqs/1.faa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | tail -n 1 | sed 's/^diamond version //')
    END_VERSIONS
    """
}
