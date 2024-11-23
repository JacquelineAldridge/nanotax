process PLOT_TAXONOMY {
    label 'process_single'
    conda "conda-forge::numpy conda-forge::seaborn conda-forge::matplotlib conda-forge::pandas"
     container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
         'ghcr.io/jacquelinealdridge/python-plots:latest' :
         'ghcr.io/jacquelinealdridge/python-plots:latest|' }"

    input:
    path summary_csv

    output:
    path "*", emit: plots
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    plot_taxonomy.py --dir_csv ${summary_csv}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
        seaborn: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('seaborn').version)")
        matplotlib: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('matplotlib').version)")
        numpy: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('numpy').version)")

    END_VERSIONS
    """

  
}
