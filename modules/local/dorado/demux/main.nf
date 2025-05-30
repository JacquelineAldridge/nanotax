process DORADO_DEMUX {
    container 'ghcr.io/dialvarezs/containers/dorado:1.0.0'

    input:
    tuple val(meta), path(basecalled_reads)

    output:
    path 'demultiplexed/*barcode*', emit: classified
    path 'demultiplexed/unclassified*', emit: unclassified
    path 'versions.yml', emit: versions

    script:
    """
    dorado demux \\
        --output-dir demultiplexed/ \\
        --emit-fastq \\
        --threads ${task.cpus} \\
        --kit-name ${params.dorado_barcoding_kit} \\
        ${basecalled_reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( dorado --version 2>&1 | tr -d '\n' )
    END_VERSIONS
    """
}
