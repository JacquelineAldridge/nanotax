process MERGE_PICRUST_OUT {
    label 'process_single'
    conda "conda-forge::pandas=2.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
                'biocontainers/pandas:2.2.1' :
                'biocontainers/pandas:2.2.1' }"
    input:
    path picrust_out_dir
    val ch_groups

    output:
    path "*.csv", emit: csv
    path "lefseinput_*.tsv", emit: lefse_input

    script:
    """
    merge_picrust_out.py --groups "${ch_groups}"

    """

}
