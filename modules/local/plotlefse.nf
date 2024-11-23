process PLOT_LEFSE {
    label 'process_single'
    conda "conda-forge::numpy conda-forge::seaborn conda-forge::matplotlib conda-forge::pandas"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'ghcr.io/jacquelinealdridge/python-plots:latest' :
            'ghcr.io/jacquelinealdridge/python-plots:latest' }"

    input:
    path lefse_results

    output:
    path "*", emit: pdf,  optional:true


    script:
    """
    plot_lefse.py -i ${lefse_results} 

    """

}
