process DORADO_BASECALLER {
    label 'process_long'

    container 'ghcr.io/dialvarezs/containers/dorado:1.0.0'
    clusterOptions "--gres=gpu:${params.dorado_gpus}"
    containerOptions { workflow.containerEngine == 'singularity' ? '--nv' : '' }

    input:
    path pod5_dir

    output:
    tuple val { [id: 'basecalled'] }, path('basecalled.ubam'), emit: reads
    path 'sequencing_summary.txt', emit: sequencing_summary
    path 'versions.yml', emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    dorado basecaller \\
        --recursive \\
        --device 'cuda:all' \\
        --trim adapters \\
        ${args} \\
        ${params.dorado_basecalling_model} \\
        ${pod5_dir} \\
    > basecalled.ubam

    dorado summary basecalled.ubam > sequencing_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( dorado --version 2>&1 | tr -d '\n' )
    END_VERSIONS
    """
}
