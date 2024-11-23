process MMSEQS_CREATE16SDB {
    label 'process_single'
    conda "bioconda::mmseqs2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:15.6f452--pl5321h6a68c12_0':
        'biocontainers/mmseqs2:15.6f452--pl5321h6a68c12_0' }"

    input:
    path(genbank_db)
    val(db_name)

    output:
    tuple val("db"),path("16S_DB/"), emit: path_db
    path "versions.yml"     , emit: versions

    script:
    if(db_name ==  'silva' )
    """
    mkdir 16S_DB
    cd 16S_DB
    mmseqs databases SILVA db tmp 
    cd ..
    """
    else if(db_name == 'genbank') 
    """
    mkdir 16S_DB
    cd 16S_DB
    mmseqs createdb ../db/16S.fna db
    mmseqs createtaxdb db tmp --ncbi-tax-dump ../db/taxonomy/ --tax-mapping-file ../db/16s.fna.taxidmapping --threads ${task.cpus}
    cd ..
        cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """
    else 
        error "invalid database, options are 'silva' or 'genbank'"

}
