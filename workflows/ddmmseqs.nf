/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { samplesheetToList } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { DIAMOND_MAKEDB_FROM_PIPE } from '../modules/local/diamond_makedb_from_pipe/main'
include { DIAMOND_BLASTP_TO_CLUSTER } from '../modules/local/diamond_blastp_to_cluster/main'
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

    DIAMOND_MAKEDB_FROM_PIPE(filelist_ch)
    diamond_db_ch = DIAMOND_MAKEDB_FROM_PIPE.out.db
    DIAMOND_BLASTP_TO_CLUSTER(
        filelist_ch.join(diamond_db_ch), 
        true, 
        params.target_cluster_size,
    )

    emit:
    versions = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
