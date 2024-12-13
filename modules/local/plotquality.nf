process PLOT_QUALITY{
    label 'process_single'
    memory { 10.GB * task.attempt }
    conda "conda-forge::numpy conda-forge::seaborn conda-forge::matplotlib conda-forge::biopython conda-forge::pandas"
     container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
         'ghcr.io/jacquelinealdridge/quality-plot-by-plots:latest' :
         'ghcr.io/jacquelinealdridge/quality-plot-by-plots:latest' }"

    input:
    file(filtered_fq)

    output:
    path "quality_plot.pdf", emit: plot

    script:
    """
    zcat ${filtered_fq} > all.fq
    length_quality_plot.py -i all.fq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}