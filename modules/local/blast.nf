process BLAST {
    label 'process_single'

    conda "bioconda::blast=2.13.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.16.0--hc155240_2':
        'biocontainers/blast:2.16.0--hc155240_2' }"

    output:
    path("db/")             , emit: db_files    
    //path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir db db/taxdb
    wget https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz && tar -xzvf 16S_ribosomal_RNA.tar.gz -C db
    wget https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz && tar -xzvf taxdb.tar.gz -C db/taxdb

    cd db

    mkdir taxonomy
    cd taxonomy
    wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/new_taxdump_2024-11-01.zip
    unzip new_taxdump_2024-11-01.zip
    cd ..
    export BLASTDB=./
    blastdbcmd -db 16S_ribosomal_RNA -entry all > 16S.fna
    blastdbcmd -db 16S_ribosomal_RNA -entry all -outfmt "%a %T" > 16s.fna.taxidmapping
    cd ..
    """
}
