/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                 } from '../modules/nf-core/multiqc/main'
include { MMSEQS_EASYSEARCH       } from '../modules/nf-core/mmseqs/easysearch'
include { SUMMARY_MMSEQS          } from '../modules/local/summarymmseqs'
include { MERGE_AND_GROUP_SAMPLES } from '../modules/local/mergeandgroupsamples'
include { PLOT_CORE               } from '../modules/local/plotcore'
include { PLOT_TAXONOMY           } from '../modules/local/plottaxonomy'
include { DIVERSITY               } from '../modules/local/diversity'
include { SEQKIT                  } from '../modules/local/seqkit'
include { SEQKIT_GREP             } from '../modules/nf-core/seqkit/grep'
include { SEQKIT_SEQ              } from '../modules/nf-core/seqkit/seq'
include { SEQKIT_FQ2FA            } from '../modules/nf-core/seqkit/fq2fa'
include { OBTAIN_IDS              } from '../modules/local/obtainids'
include { PICRUST2                } from '../modules/local/picrust2'
include { MERGE_PICRUST_OUT       } from '../modules/local/mergepicrustout'
include { LEFSE                   } from '../modules/local/lefse'
include { PLOT_LEFSE              } from '../modules/local/plotlefse'

include { BASECALLING             } from '../subworkflows/local/basecalling/main'
include { PREPARE_DATABASES       } from '../subworkflows/local/prepare_databases/main'
include { QUALITY_CONTROL         } from '../subworkflows/local/quality_control/main'
include { TAXONOMIC_ASSIGNMENT    } from '../subworkflows/local/taxonomic_assignment/main'

include { paramsSummaryMap        } from 'plugin/nf-schema'
include { paramsSummaryMultiqc    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText  } from '../subworkflows/local/utils_nfcore_nanotax_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NANOTAX {
    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    /*
     * Database preparation
     */
    PREPARE_DATABASES(
        params.emu_database,
        params.skip_emu,
        params.skip_mmseqs2,
    )

    /*
     * Basecalling and demultiplexing using Dorado
     */
    if (!params.skip_basecalling) {
        ch_pod5_dir = channel.value(
            [[id: 'reads'], file(params.dorado_pod5_dir, checkIfExists: true, type: 'dir')]
        )

        BASECALLING(ch_pod5_dir, ch_samplesheet, params.dorado_barcoding_kit)

        ch_samples = BASECALLING.out.samples
        ch_versions = ch_versions.mix(BASECALLING.out.versions)
    }
    else {
        ch_samples = ch_samplesheet
    }

    /*
     * Quality control and filtering
     */
    if (!params.skip_qc) {
        QUALITY_CONTROL(ch_samples, params.filtlong_sampling)
    }

    /*
     * Taxonomic assignment
     */
    TAXONOMIC_ASSIGNMENT(
        QUALITY_CONTROL.out.reads,
        PREPARE_DATABASES.out.emu_database,
        PREPARE_DATABASES.out.mmseqs2_database,
        params.skip_emu,
        params.skip_mmseqs2,
    )

    // // Taxonomic assignment
    // if(params.mmseqs2_download_db && params.mmseqs2_db_name == 'genbank'){
    //     BLASTCMD()
    //     MMSEQS_CREATE16SDB(BLASTCMD.out.db_files,params.mmseqs2_db_name)

    // }else if(params.mmseqs2_download_db && params.mmseqs2_db_name == 'silva'){
    //     MMSEQS_CREATE16SDB([],params.mmseqs2_db_name)
    // Plots for Taxonomic assignment
    // ch_groups = ch_samplesheet.map{meta,path-> "${meta.id}:${meta.group}"}.collect()
    // PLOT_TAXONOMY((MERGE_AND_GROUP_SAMPLES.out.csv_sample.mix(MERGE_AND_GROUP_SAMPLES.out.csv_group)).flatten()) //csv_group
    // PLOT_CORE(MERGE_AND_GROUP_SAMPLES.out.csv_core,SUMMARY_MMSEQS.out.taxlineage.collect(),ch_groups)

    // // Diversity
    // // ToDo: Solo si hay grupos
    // if(!params.skip_diversity){
    //     ch_groups_info_all = MMSEQS_EASYSEARCH.out.tsv.map{meta,tsv -> "${meta.id}:${meta.group}:${meta.subgroup}:${meta.subsubgroup}"}.collect()
    //     DIVERSITY(MERGE_AND_GROUP_SAMPLES.out.csv_div_nreads,ch_groups_info_all)//ch_groups)
    //     ch_versions = ch_versions.mix(DIVERSITY.out.versions.first())
    // }
    // Functional prediction
    // if(params.skip_functional_prediction){

    // }else if(!params.mmseqs2_download_db){
    //     print("ToDo: completar")
    //     ch_db_dir = Channel.fromPath(params.mmseqs2_db_dir)
    // }
    // MMSEQS_EASYSEARCH(ch_input_tax,MMSEQS_CREATE16SDB.out.path_db)
    // ch_versions = ch_versions.mix(MMSEQS_EASYSEARCH.out.versions.first())

    // ch_mmseqs_output = MMSEQS_EASYSEARCH.out.tsv//.map{meta,tsv -> tsv}.collect()
    // //ch_first_group = MMSEQS_EASYSEARCH.out.tsv.map{meta,tsv -> meta}.first()
    // ch_groups_info = MMSEQS_EASYSEARCH.out.tsv.map{meta,tsv -> "${meta.id}:${meta.group}"}.collect()
    // ch_samples = ch_samplesheet.map{meta,path-> "${meta.id}"}.collect()
    // SUMMARY_MMSEQS(ch_mmseqs_output,ch_samples)
    // MERGE_AND_GROUP_SAMPLES(SUMMARY_MMSEQS.out.summary_csv.collect())//, SUMMARY_MMSEQS.out.abundance_picrust.collect())

    // // Plots for Taxonomic assignment
    // ch_groups = ch_samplesheet.map{meta,path-> "${meta.id}:${meta.group}"}.collect()
    // PLOT_TAXONOMY((MERGE_AND_GROUP_SAMPLES.out.csv_sample.mix(MERGE_AND_GROUP_SAMPLES.out.csv_group)).flatten()) //csv_group
    // PLOT_CORE(MERGE_AND_GROUP_SAMPLES.out.csv_core,SUMMARY_MMSEQS.out.taxlineage.collect(),ch_groups)

    // // Diversity
    // // ToDo: Solo si hay grupos
    // if(!params.skip_diversity){
    //     ch_groups_info_all = MMSEQS_EASYSEARCH.out.tsv.map{meta,tsv -> "${meta.id}:${meta.group}:${meta.subgroup}:${meta.subsubgroup}"}.collect()
    //     DIVERSITY(MERGE_AND_GROUP_SAMPLES.out.csv_div_nreads,ch_groups_info_all)//ch_groups)
    //     ch_versions = ch_versions.mix(DIVERSITY.out.versions.first())
    // }
    // // Functional prediction
    // if(params.skip_functional_prediction){
    //     ch_input_picrust = (SUMMARY_MMSEQS.out.abundance_picrust.join(ch_input_tax)).map{meta,tsv,fastq -> [tsv,fastq]}
    //     SEQKIT(ch_input_picrust) //ch_input_tax.map{meta, path -> path}.collect(),SUMMARY_MMSEQS.out.abundance_picrust.collect())
    //     (SEQKIT.out.abundance, SEQKIT.out.fasta)
    //     ch_versions = ch_versions.mix(PICRUST2.out.versions.first())

    //     MERGE_PICRUST_OUT(PICRUST2.out.dir.collect(), ch_groups)
    //     LEFSE(MERGE_PICRUST_OUT.out.lefse_input.flatten())
    //     PLOT_LEFSE(LEFSE.out.lefse_output.flatten())
    //     // ToDo:LEFSE SOLO SI HAY GRUPOS
    // }
    //     OBTAIN_IDS(SUMMARY_MMSEQS.out.abundance_picrust)
    //     ch_input_seqkit =  OBTAIN_IDS.out.abundance.join(ch_input_tax)
    //                         .join(OBTAIN_IDS.out.ids)
    //                         .multiMap{meta,tsv,fastq,ids ->
    //                             sequences: [meta, fastq]
    //                             ids: [ids]
    //     }
    //     SEQKIT_GREP(ch_input_seqkit.sequences,ch_input_seqkit.ids)
    //     SEQKIT_FQ2FA(SEQKIT_GREP.out.filter)
    //     SEQKIT_SEQ(SEQKIT_FQ2FA.out.fasta)
    //     ch_input_picrust = OBTAIN_IDS.out.abundance.join( SEQKIT_SEQ.out.fastx)
    //                         .multiMap{meta,abundance_table,fasta->
    //                         abundance_table: [meta, abundance_table]
    //                         fasta: [fasta]
    //                         }
    //     PICRUST2(ch_input_picrust.abundance_table,ch_input_picrust.fasta)
    //     ch_versions = ch_versions.mix(PICRUST2.out.versions.first())
    //     MERGE_PICRUST_OUT(PICRUST2.out.dir.collect(), ch_groups)
    //     LEFSE(MERGE_PICRUST_OUT.out.lefse_input.flatten())

    //     // ToDo:LEFSE: SOLO SI HAY GRUPOS ; version
    // }

    // Collate and save software versions

    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nanotax_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )
    ch_multiqc_custom_config = params.multiqc_config
        ? channel.fromPath(params.multiqc_config, checkIfExists: true)
        : channel.empty()
    ch_multiqc_logo = params.multiqc_logo
        ? channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : channel.empty()

    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )
    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true,
        )
    )

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
