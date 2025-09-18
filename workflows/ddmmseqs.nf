/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { samplesheetToList } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { DIAMOND_MAKEDB_FROM_PIPE } from '../modules/local/diamond_makedb_from_pipe/main'
include { DIAMOND_BLASTP_TO_CLUSTER } from '../modules/local/diamond_blastp_to_cluster/main'
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

    filelist_ch = manifest.map {
        name, filelist ->
        [
            ['id': name],
            file(filelist),
        ]
    }
    
    // Greedy clustering to generate chunks
    DIAMOND_MAKEDB_FROM_PIPE(filelist_ch)
    diamond_db_ch = DIAMOND_MAKEDB_FROM_PIPE.out.db
    DIAMOND_BLASTP_TO_CLUSTER(
        filelist_ch.join(diamond_db_ch), 
        true, 
        params.target_cluster_size,
    )
    def idx = 0
    seq_chunks_ch = DIAMOND_BLASTP_TO_CLUSTER.out.cluster_seqs
        .map{ meta, seqs -> 
            [meta + ['n_seqs': seqs.size()], seqs] 
        }
        .transpose()
        .map{
            meta, seq ->
            idx += 1
            [meta + ['idx': idx], seq]
        }

    // clustering of chunks
    MMSEQS_CREATEDB(seq_chunks_ch)
    MMSEQS_LINCLUST(MMSEQS_CREATEDB.out.db)
    joined_db_clustered_db_ch = MMSEQS_CREATEDB.out.db
        .join(MMSEQS_LINCLUST.out.db_cluster)
        .multiMap{ meta, db, db_cluster -> 
            db: [meta, db]
            db_cluster: [meta, db_cluster]
        }
    MMSEQS_CREATETSV(
        joined_db_clustered_db_ch.db_cluster,
        joined_db_clustered_db_ch.db,
        joined_db_clustered_db_ch.db,
    )

    // concatentate results
    tsv_concat_ch = MMSEQS_CREATETSV.out.tsv
        .map{ meta, tsv -> [groupKey(meta.remove('idx'), meta.n_seqs), tsv] }
        .groupTuple()
        .map{ meta, tsvs ->
            [meta, meta.id, tsvs]}
    CONCATENATE(tsv_concat_ch)

    emit:
    clusters = CONCATENATE.out.concatenated_file
    versions = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
