process DORADO_DEMUX {
    label 'process_medium'

    container 'ghcr.io/dialvarezs/containers/dorado:1.0.0'

    input:
    tuple val(meta), path(basecalled_reads)
    path sample_sheet
    val kit_name

    output:
    tuple val(meta), path('demultiplexed/*[!unclassified].fastq'), emit: classified
    tuple val(meta), path('demultiplexed/*unclassified.fastq'), emit: unclassified
    tuple val(meta), path('demultiplexed/barcoding_summary.txt'), emit: summary
    path 'versions.yml', emit: versions

    script:
    def args = task.ext.args ?: ''
    def sample_sheet_arg = sample_sheet ? "--sample-sheet ${sample_sheet}" : ""
    """
    dorado demux \\
        --output-dir demultiplexed \\
        --emit-fastq \\
        --emit-summary \\
        --threads ${task.cpus} \\
        --kit-name ${kit_name} \\
        ${sample_sheet_arg} \\
        ${basecalled_reads} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$( dorado --version 2>&1 )
    END_VERSIONS
    """
}
