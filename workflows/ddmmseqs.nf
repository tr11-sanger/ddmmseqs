/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { samplesheetToList } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { DIAMOND_MAKEDB_FROM_PIPE } from '../modules/local/diamond_makedb_from_pipe/main'
include { DIAMOND_BLASTP_TO_MST } from '../modules/local/diamond_blastp_to_mst/main'
include { MST_TO_CLUSTER } from '../modules/local/mst_to_cluster/main'
include { MMSEQS_CREATEDB } from '../modules/nf-core/mmseqs/createdb/main'
include { MMSEQS_LINCLUST } from '../modules/nf-core/mmseqs/linclust/main'
include { MMSEQS_CREATETSV } from '../modules/nf-core/mmseqs/createtsv/main'
include { CONCATENATE } from '../modules/local/concatenate/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DDMMSEQS {
    main:
    ch_versions = Channel.empty()

    // Parse samplesheet and fetch reads
    manifest = Channel.fromList(samplesheetToList(params.samplesheet, "${workflow.projectDir}/assets/schema_input.json"))
    def manifest_idx = 0
    faa_chunks_ch = manifest.map {
        collection_name, chunk_name, filelist ->
        manifest_idx += 1
        return [
            ['id': "${collection_name}_${chunk_name}", 'idx': manifest_idx, 
             'collection_id': collection_name, 'chunk_id': chunk_name],
            file(filelist),
        ]
    }
    
    // Insert number of chunks in each collection
    faa_chunks_ch = faa_chunks_ch
        .map { meta, filelist -> [['id': meta.collection_id], [meta, filelist]]}
        .groupTuple()
        .map { k, vs ->
            def n = vs.size() 
            return [k, vs.collect{meta, filelist -> [meta + ['n_chunks': n], filelist]}] }
        .transpose()
        .map { _k, v -> 
            def (meta, filelist) = v 
            return [meta, filelist]
        }
    faa_chunks_ch.view{ "faa_chunks_ch - ${it}" }
    
    // Greedy clustering to generate chunks
    DIAMOND_MAKEDB_FROM_PIPE(faa_chunks_ch)
    diamond_db_ch = DIAMOND_MAKEDB_FROM_PIPE.out.db

    diamond_blastp_ch = faa_chunks_ch
        .combine(diamond_db_ch)
        .filter{ meta1, _filelist, meta2, _db -> 
            (meta1.collection_id==meta2.collection_id) && (meta1.idx<=meta2.idx)
        }
        .map{ meta1, filelist, meta2, db -> 
            def meta = ['id': "${meta1.collection_id}_${meta1.chunk_id}_${meta2.chunk_id}", 
                        'collection_id': meta1.collection_id, 'n_chunks': meta1.n_chunks]
            return [['id': meta.collection_id], [meta, filelist, db]] 
        } 
        .groupTuple()
        .map { k, vs ->
            def n = vs.size() 
            return [k, vs.collect{meta, filelist, db -> [meta + ['n_perm': n], filelist, db]}] }
        .transpose()
        .map { k, v -> v }
    diamond_blastp_ch.view{ "diamond_blastp_ch - ${it}" }

    DIAMOND_BLASTP_TO_MST(
        diamond_blastp_ch,
        false, 
    )
    cluster_ch = DIAMOND_BLASTP_TO_MST.out.nodes
        .join( DIAMOND_BLASTP_TO_MST.out.linkage )
        .join( diamond_blastp_ch.map { meta, filelist, _db -> [meta, filelist] } )
        .map { meta, nodes, linkage, filelist -> [groupKey(['id': meta.collection_id], meta.n_perm), [nodes, linkage, filelist]] }
        .groupTuple()
        .map { meta, vs -> 
            def (nodes, linkages, filelists) = vs.transpose()
            return [meta, [nodes, linkages].transpose(), filelists]
        }
    cluster_ch.view{ "cluster_ch - ${it}" }
    MST_TO_CLUSTER(cluster_ch, params.target_cluster_size)

    def idx = 0
    seq_chunks_ch = MST_TO_CLUSTER.out.cluster_seqs
        .map { meta, seqs -> 
            [meta + ['n_seqs': seqs.size()], seqs] 
        }
        .transpose()
        .map{
            meta, seq ->
            idx += 1
            [meta + ['idx': idx], seq]
        }
    seq_chunks_ch.view{ "seq_chunks_ch - ${it}" }

    // clustering of chunks
    MMSEQS_CREATEDB(seq_chunks_ch)
    MMSEQS_LINCLUST(MMSEQS_CREATEDB.out.db)
    joined_db_clustered_db_ch = MMSEQS_CREATEDB.out.db
        .join(MMSEQS_LINCLUST.out.db_cluster)
        .multiMap { meta, db, db_cluster -> 
            db: [meta, db]
            db_cluster: [meta, db_cluster]
        }
    joined_db_clustered_db_ch.db_cluster.view{ "joined_db_clustered_db_ch.db_cluster - ${it}" }
    MMSEQS_CREATETSV(
        joined_db_clustered_db_ch.db_cluster,
        joined_db_clustered_db_ch.db,
        joined_db_clustered_db_ch.db,
    )

    // concatentate results
    tsv_concat_ch = MMSEQS_CREATETSV.out.tsv
        .map { meta, tsv -> 
            [groupKey(['id': meta.id], meta.n_seqs), tsv] 
        }
        .groupTuple()
        .map { meta, tsvs -> [meta, meta.id, tsvs] }
    tsv_concat_ch.view{ "tsv_concat_ch - ${it}" }
    CONCATENATE(tsv_concat_ch)

    // output representative sequences?

    emit:
    clusters = CONCATENATE.out.concatenated_file
    versions = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
