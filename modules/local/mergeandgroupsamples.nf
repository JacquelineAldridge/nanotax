
process MERGE_AND_GROUP_SAMPLES {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'ghcr.io/dialvarezs/containers/polars:1.3.0' :
            'ghcr.io/dialvarezs/containers/polars:1.3.0' }"

    input:

    path summary_by_sample

    output:
    path "group/", emit: csv_group,optional: true
    path "sample/", emit: csv_sample
    path "sample/species.csv", emit: csv_species
    //path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    mkdir group
    mkdir sample
    merge_taxonomies.py --csv "${summary_by_sample}" --db ${params.taxonomic_assignament.db_name}

    """


}
