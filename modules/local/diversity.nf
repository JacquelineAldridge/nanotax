process DIVERSITY {
    label 'process_single'
    conda "conda-forge::r-base conda-forge::r-tidyverse conda-forge::r-vegan  conda-forge::r-ggplot2 conda-forge::r-cowplot conda-forge::r-optparse"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'ghcr.io/jacquelinealdridge/diversity:latest' :
            'ghcr.io/jacquelinealdridge/diversity:latest' }"

    input:
    path(species_nreads)
    val(groups)

    output:
    path "*.csv", emit: diversity_csv
    path "*.pdf", emit: plot
    path "versions.yml", emit: versions

    script:   
    """
    diversity.R ${species_nreads} "${groups}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n 1 | awk '{print \$3}' )
    END_VERSIONS
    """

}
