
process LEFSE {
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/lefse:1.1.2--pyhdfd78af_0':
        'biocontainers/lefse:1.1.2--pyhdfd78af_0' }"

    input:
    path input_file

    output:
    path "${name}/", emit: lefse_output



    script: 
    name = input_file.simpleName - ~/lefseinput_/
    """
    mkdir ${name}
    lefse_format_input.py ${input_file} ${name}_format.in -c 1 -s -1 -u 2 -o 1000000
    lefse_run.py ${name}_format.in  ${name}_lefse_results.csv 
    lefse_plot_res.py ${name}_lefse_results.csv  ${name}.res.png
    cp ${name}_lefse_results.csv ${name}
    """
}
