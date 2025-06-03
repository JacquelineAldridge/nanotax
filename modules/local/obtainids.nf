process OBTAIN_IDS {
    label 'process_single'
    debug true
    input:
    tuple val(meta), path(abundance_table)

    output:
    tuple val(meta), path (abundance_table), emit: abundance
    tuple val(meta), path("${meta.id}_ids.txt"), emit: ids
    // path "versions.yml"           , emit: versions

    script:    
    """
    cat ${abundance_table} | tail -n +2 | awk -F '\t' '{print \$1}' > ${meta.id}_ids.txt
    """
    // touch versions.yml

    // stub:
    // def args = task.ext.args ?: ''
    // def prefix = task.ext.prefix ?: "${meta.id}"
    // // fasta or fastq. Exact pattern match .fasta or .fa suffix with optional .gz (gzip) suffix
    // // def suffix = task.ext.suffix ?: "${sequence}" ==~ /(.*f[astn]*a(.gz)?$)/ ? "fa" : "fq"
    // """
    // echo "" > ids.txt
    // """
        // touch versions.yml

}
