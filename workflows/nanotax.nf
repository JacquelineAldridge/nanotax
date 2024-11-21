/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { NANOQ as NANOQ_FILTER  } from '../modules/nf-core/nanoq/main'
include { NANOQ as NANOQ_QC_RAW  } from '../modules/nf-core/nanoq/main'
include { FILTLONG               } from '../modules/nf-core/filtlong/main'
include { BLAST as BLASTCMD      } from '../modules/local/blast'
include { MMSEQS_CREATE16SDB     } from '../modules/local/mmseqs/create16sdb'
include { MMSEQS_EASYSEARCH      } from '../modules/nf-core/mmseqs/easysearch'
include { SUMMARY_MMSEQS         } from '../modules/local/summarymmseqs'
include { MERGE_AND_GROUP_SAMPLES } from '../modules/local/mergeandgroupsamples'
include { PLOT_CORE              } from '../modules/local/plotcore'
include { PLOT_TAXONOMY          } from '../modules/local/plottaxonomy'

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_nanotax_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NANOTAX {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    
    // BASECALLING AND DEMUX
    if(params.qc.basecalling){
        // ToDo: cambiar por basecalling
        ch_input_qc = ch_samplesheet
    }else{
        ch_input_qc = ch_samplesheet
    }
    

    // QC
    if(params.qc.run){
        FASTQC (
            ch_input_qc
        )

        NANOQ_FILTER(ch_input_qc,'fastq.gz')
        NANOQ_QC_RAW(ch_input_qc,'fastq.gz')

        if(params.qc.subsampling>0){
            ch_mix = NANOQ_FILTER.out.reads.map{it -> [it[0],[],it[1]]}
            FILTLONG(ch_mix)
            ch_input_tax = FILTLONG.out.reads
            ch_versions.mix(FILTLONG.out.versions.first())
        }else{
            ch_input_tax = NANOQ_FILTER.out.reads

        }

        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]},NANOQ_QC_RAW.out.stats.collect{it[1]},NANOQ_FILTER.out.stats.collect{it[1]})
        ch_versions = ch_versions.mix(FASTQC.out.versions.first(),NANOQ_FILTER.out.versions.first())
    
    }else{
        ch_input_tax = ch_input_qc
    }

    // Taxonomic assignment
    if(params.taxonomic_assignament.download_db && params.taxonomic_assignament.db_name == 'genbank'){
        BLASTCMD()
        MMSEQS_CREATE16SDB(BLASTCMD.out.db_files,params.taxonomic_assignament.db_name)

    }else if(params.taxonomic_assignament.download_db && params.taxonomic_assignament.db_name == 'silva'){
        MMSEQS_CREATE16SDB([],params.taxonomic_assignament.db_name)
    }else if(!params.taxonomic_assignament.download_db){
        print("ToDo: completar")
        ch_db_dir = Channel.fromPath(params.taxonomic_assignament.db_dir)
    }
    MMSEQS_EASYSEARCH(ch_input_tax,MMSEQS_CREATE16SDB.out.path_db)
    ch_versions = ch_versions.mix(MMSEQS_EASYSEARCH.out.versions.first())

    ch_mmseqs_output = MMSEQS_EASYSEARCH.out.tsv//.map{meta,tsv -> tsv}.collect()
    //ch_first_group = MMSEQS_EASYSEARCH.out.tsv.map{meta,tsv -> meta}.first()
    ch_groups_info = MMSEQS_EASYSEARCH.out.tsv.map{meta,tsv -> "${meta.id}:${meta.group}"}.collect()
    SUMMARY_MMSEQS(ch_mmseqs_output)
    MERGE_AND_GROUP_SAMPLES(SUMMARY_MMSEQS.out.summary_csv.collect())

    // Plots for Taxonomic assignment
    PLOT_TAXONOMY((MERGE_AND_GROUP_SAMPLES.out.csv_sample.mix(MERGE_AND_GROUP_SAMPLES.out.csv_group)).flatten()) //csv_group
    
    ch_groups = ch_samplesheet.map{meta,path-> "${meta.id}:${meta.group}"}.collect()


    PLOT_CORE(MERGE_AND_GROUP_SAMPLES.out.csv_species,SUMMARY_MMSEQS.out.taxlineage.collect(),ch_groups)
    // Diversity
    // Functional prediction
    // Write on db
    


    // Collate and save software versions
    
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  ''  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
