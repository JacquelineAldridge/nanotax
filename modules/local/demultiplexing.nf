process DEMULTIPLEXING {
    label 'dorado'

    input:
      tuple val(name), path(basecalled_reads)

    output:
      path('demultiplexed/*barcode*')    , emit: classified
      path('demultiplexed/unclassified*'), emit: unclassified

    script:    
    """
    dorado demux \
      --output-dir demultiplexed/ \
      --emit-fastq \
      --threads ${task.cpus} \
      --kit-name ${params.basecalling.barcoding_kit} \
      ${basecalled_reads}

    """

}
