
process MERGE_AND_GROUP_SAMPLES {
    label 'process_single'
    conda "conda-forge::polars=1.14.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'ghcr.io/dialvarezs/containers/polars:1.3.0' :
            'ghcr.io/dialvarezs/containers/polars:1.3.0' }"

    input:

    path summary_by_sample

    output:
    path "group/", emit: csv_group,optional: true
    path "sample/", emit: csv_sample
    path "core/last_assignment.csv", emit: csv_core
    path "diversity_nreads.csv", emit: csv_div_nreads

    //path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    mkdir group core sample
    merge_taxonomies.py --csv "${summary_by_sample}" --db ${params.taxonomic_assignment.db_name}

    """


}
