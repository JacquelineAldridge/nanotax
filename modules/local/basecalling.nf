process BASECALLING {
    label 'dorado'
    clusterOptions "--gres=gpu:${params.dorado_gpus}"
    cpus { 4 * params.dorado_gpus }
    memory "${16 * params.dorado_gpus}G"

    input:
      path(pod5_dir)

    output:
      tuple val('basecalled'), path('basecalled.ubam'), emit: reads
      path('sequencing_summary.txt')                  , emit: sequencing_summary
      //path "versions.yml"           , emit: versions


    script:
    def args = task.ext.args ?: ''
    """
    dorado basecaller \
    --recursive \
    --device 'cuda:all' \
    ${args} \
    ${params.dorado_basecalling_model} \
    ${pod5_dir} \
  > basecalled.ubam

  dorado summary basecalled.ubam > sequencing_summary.txt
    """
}
