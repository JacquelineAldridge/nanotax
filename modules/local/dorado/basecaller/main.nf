process DORADO_BASECALLER {
    label 'process_long'
    label 'process_gpu'

    container 'ghcr.io/dialvarezs/containers/dorado:1.0.0'
    clusterOptions "--gres=gpu:${params.dorado_basecalling_gpus}"

    input:
    tuple val(meta), path(pod5_dir)

    output:
    tuple val(meta), path("${prefix}_basecalled.ubam"), emit: reads
    tuple val(meta), path("${prefix}_sequencing_summary.txt"), emit: sequencing_summary
    path 'versions.yml', emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    dorado basecaller \\
        --recursive \\
        --device 'cuda:all' \\
        --trim adapters \\
        ${args} \\
        ${params.dorado_basecalling_model} \\
        ${pod5_dir} \\
    > ${prefix}_basecalled.ubam

    dorado summary ${prefix}_basecalled.ubam > ${prefix}_sequencing_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$( dorado --version 2>&1 )
    END_VERSIONS
    """
}
