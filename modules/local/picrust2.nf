
process PICRUST2 {
    label 'process_high'
    conda "bioconda::picrust2=2.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'biocontainers/picrust2:2.5.2--pyhdfd78af_0' :
        'biocontainers/picrust2:2.5.2--pyhdfd78af_0' }"

    input:
     file abundances_tsv 
     file fasta

    output:
     path("${name}"), emit: dir
     path "versions.yml", emit: versions
     
    script:
    name = abundances_tsv.simpleName - ~/reads_/
    """
    picrust2_pipeline.py -s ${fasta} -i ${abundances_tsv} -o ${name} --stratified -p ${task.cpus}
    rm -r ${name}/intermediate/
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picrust2: \$(picrust2_pipeline.py --version | sed 's/picrust2_pipeline.py //g')
    END_VERSIONS
    """


    
}
