process SUMMARY_MMSEQS {
    memory { 6.GB * task.attempt }
    label 'process_single'
    conda "conda-forge::polars=1.14.0"
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 3
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'ghcr.io/dialvarezs/containers/polars:1.3.0' :
            'ghcr.io/dialvarezs/containers/polars:1.3.0' }"

    input:
    tuple val(meta), path(mmseqs_tsv)
    val(samples)
    //val(meta_groups)

    output:
    path("*.csv"), emit: summary_csv
    path("taxlineage/${meta.id}_taxlineage.csv"), emit: taxlineage
    tuple val(meta), path("reads_*.tsv"), emit: abundance_picrust

    //ToDo: polars version
    //path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    mkdir taxlineage
    summary_mmseqs.py --db ${params.taxonomic_assignment.db_name} --mmseqs_tsv ${mmseqs_tsv} --min_aln ${params.taxonomic_assignment.min_aln} --min_identity ${params.taxonomic_assignment.min_identity} --group ${meta.group} --sample ${meta.id} 

    """
/*
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        summarymmseqs: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        summarymmseqs: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
    */
}
