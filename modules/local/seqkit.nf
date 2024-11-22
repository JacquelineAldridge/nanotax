process SEQKIT {
    label 'process_single'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.8.1--h9ee0642_0':
        'biocontainers/seqkit:2.8.1--h9ee0642_0' }"

    input:
    tuple path(abundance_table), path(fastq) 
    //path bam

    output:
    path (abundance_table), emit: abundance
    path "filtered_i.fasta", emit: fasta
    path "versions.yml"           , emit: versions

    script:    
    """
    cat ${abundance_table} | tail -n +2 | awk -F '\t' '{print \$1}' > ids.txt
    seqkit grep -f ids.txt ${fastq} > filtered.fastq
    seqkit fq2fa filtered.fastq -o filtered.fasta
    seqkit seq -i filtered.fasta > filtered_i.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // fasta or fastq. Exact pattern match .fasta or .fa suffix with optional .gz (gzip) suffix
    def suffix = task.ext.suffix ?: "${sequence}" ==~ /(.*f[astn]*a(.gz)?$)/ ? "fa" : "fq"

    """
    echo "" | gzip > ${prefix}.${suffix}.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
    END_VERSIONS
    """
}
