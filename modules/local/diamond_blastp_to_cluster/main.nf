process DIAMOND_BLASTP_TO_CLUSTER {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    //     ? 'https://depot.galaxyproject.org/singularity/diamond:2.1.12--hdb4b4cc_1'
    //     : 'biocontainers/diamond:2.1.12--hdb4b4cc_1'}"

    input:
    tuple val(meta), path(filelist), path(db)
    val save_paf

    output:
    tuple val(meta), path('*_nodes.{tsv,tsv.gz}'), emit: nodes
    tuple val(meta), path('*_linkage.{tsv,tsv.gz}'), emit: linkage
    tuple val(meta), path('*_clusters.{tsv,tsv.gz}'), emit: clusters
    tuple val(meta), path('*_cluster_seqs/*.{faa,faa.gz}'), emit: cluster_seqs
    tuple val(meta), path('*.{paf,paf.gz}'), optional: true, emit: paf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def outfmt = 103
    def blast_columns = false
    def columns = blast_columns ? "${blast_columns}" : ''
    def tee_cmd = save_paf ? "| tee ${prefix}.paf" : ''

    """
    cat ${filelist} \\
    | while read fp; do
        if [[ \$fp =~ \\.gz\$ ]]; then
            gunzip -c \$fp
        else
            cat \$fp
        fi
    done \\
    | diamond \\
        blastp \\
        --threads ${task.cpus} \\
        --db ${db} \\
        --outfmt ${outfmt} ${columns} \\
        ${args} \\
    ${tee_cmd} | python ${moduleDir}/bin/greedy_cluster.py \\
        --filelist ${filelist} \\
        --out_nodes ${prefix}_nodes.tsv.gz \\
        --out_linkage ${prefix}_linkage.tsv.gz \\
        --out_clusters ${prefix}_clusters.tsv.gz \\
        --out_cluster_seqs_dir ${prefix}_cluster_seqs
        ${args2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | tail -n 1 | sed 's/^diamond version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_nodes.txt.gz
    touch ${prefix}_linkage.tsv.gz
    touch ${prefix}_clusters.tsv.gz
    mkdir ${prefix}_cluster_seqs
    touch ${prefix}_cluster_seqs/1.faa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | tail -n 1 | sed 's/^diamond version //')
    END_VERSIONS
    """
}
