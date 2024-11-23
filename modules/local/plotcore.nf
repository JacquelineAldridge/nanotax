
process PLOT_CORE {
    label 'process_single'
    conda "conda-forge::python-kaleido=0.2.1 conda-forge::plotly=5.17.0 conda-forge::pandas=2.1.1"

    //conda "conda-forge::pandas=2.1.1 conda-forge::plotly=5.17.0 conda-forge::python-kaleido=0.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'ghcr.io/jacquelinealdridge/python-kaleido:latest' :
            'ghcr.io/jacquelinealdridge/python-kaleido:latest' }"

    input:
    path mmseqs_tsv
    path taxlineage
    val(groups)

    output:
     path "*.pdf"
     path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    plot_core.py --taxlineage_csv '${taxlineage}' --mmseqs_tsv ${mmseqs_tsv} --groups '${groups}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')



    END_VERSIONS
    """
}
