process PIGZ {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/dialvarezs/containers/utils:latest' :
        'ghcr.io/dialvarezs/containers/utils:latest' }"


    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq_comp

    script:
    """
    pigz -p ${task.cpus} -c ${fastq} > ${meta.id}.fastq.gz
    """

}
