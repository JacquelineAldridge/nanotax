process BASECALLING_FILTERING {
    label 'process_single'

    //conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
      tuple val(name), path(ubam)
    
    output:
      tuple val('reads_pass'), path('*.pass.ubam'), emit: reads_pass
      tuple val('reads_fail'), path('*.fail.ubam'), emit: reads_fail
      //path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    samtools view \
    -e '[qs] >= ${params.basecalling.qscore_filter}' ${ubam} \
    --output ${ubam.baseName}.pass.ubam \
    --unoutput ${ubam.baseName}.fail.ubam \
    --bam \
    --threads ${task.cpus} 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        basecallingfiltering: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

}
