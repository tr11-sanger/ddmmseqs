process MST_TO_CLUSTER {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    //     ? 'https://depot.galaxyproject.org/singularity/diamond:2.1.12--hdb4b4cc_1'
    //     : 'biocontainers/diamond:2.1.12--hdb4b4cc_1'}"

    input:
    tuple val(meta), val(node_linkage_lists), val(filelists)
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

    def node_linkage_fp = "node_linkage.txt"
    def node_linkage_echo = node_linkage_lists.collect{ n, l -> "echo ${n}\t${l} >> ${node_linkage_fp}" }.join("\n")

    def filelists_fp = "filelists.txt"
    def filelists_echo = filelists.collect{ v -> "echo ${v} >> ${filelists_fp}" }.join("\n")

    """
    :> ${node_linkage_fp}
    ${node_linkage_echo}

    :> ${filelists_fp}
    ${filelists_echo}

    python ${moduleDir}/bin/mst_to_cluster.py \\
        --linkagelists ${node_linkage_fp} \\
        --filelists ${filelists_fp} \\
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
