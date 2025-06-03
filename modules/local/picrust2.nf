
process PICRUST2 {
    label 'process_high'
    debug true
    conda "bioconda::picrust2=2.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'biocontainers/picrust2:2.5.2--pyhdfd78af_0' :
        'biocontainers/picrust2:2.5.2--pyhdfd78af_0' }"

    input:
     tuple val(meta),path (abundance_table)
     path(fasta)

    output:
     path("${name}"), emit: dir
     path "versions.yml", emit: versions
     
    script:
    name = abundance_table.simpleName - ~/reads_/
    fasta = fasta.baseName - ~/.gz/
    """
    gzip -d ${fasta}.gz
    
    picrust2_pipeline.py -s ${fasta} -i ${abundance_table} -o ${name} --stratified -p ${task.cpus}
    rm -r ${name}/intermediate/
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picrust2: \$(picrust2_pipeline.py --version | sed 's/picrust2_pipeline.py //g')
    END_VERSIONS
    """
    
}
