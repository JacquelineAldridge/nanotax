include { FASTQC                } from '../../../modules/nf-core/fastqc/main'
include { FILTLONG              } from '../../../modules/nf-core/filtlong/main'
include { NANOQ as NANOQ_FILTER } from '../../../modules/nf-core/nanoq/main'
include { NANOQ as NANOQ_QC     } from '../../../modules/nf-core/nanoq/main'
include { PLOT_QUALITY          } from '../../../modules/local/plotquality'


workflow QUALITY_CONTROL {
    take:
    ch_reads        // channel: [ val(meta), path(fastq) ]
    filtlong_sampling // float: sampling rate for filtlong, 0 to disable

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()


    NANOQ_FILTER(ch_reads, 'fastq.gz')

    ch_reads_filtered = NANOQ_FILTER.out.reads
    if (filtlong_sampling > 0) {
        FILTLONG(
            NANOQ_FILTER.out.reads.map { meta, fastq -> [meta, [], fastq] }
        )
        ch_reads_filtered = FILTLONG.out.reads

        ch_versions = ch_versions.mix(FILTLONG.out.versions.first())
    }

    // Adds '-qc' suffix to the id, so that it can be distinguished from the original reads in FastQC and nanoq.
    ch_reads_filtered_suffix = ch_reads_filtered.map { meta, fastq ->
        [meta + [id: "${meta.id}-qc"], fastq]
    }

    NANOQ_QC(ch_reads.mix(ch_reads_filtered_suffix), 'fastq.gz')
    FASTQC(ch_reads.mix(ch_reads_filtered_suffix))

    PLOT_QUALITY(NANOQ_FILTER.out.reads.collect { it[1] })


    ch_versions = ch_versions.mix(
        FASTQC.out.versions.first(),
        NANOQ_QC.out.versions.first(),
        NANOQ_FILTER.out.versions.first(),
    )
    ch_multiqc_files = ch_multiqc_files.mix(
        FASTQC.out.zip.collect { it[1] },
        NANOQ_QC.out.stats.collect { it[1] },
        NANOQ_FILTER.out.stats.collect { it[1] },
    )

    emit:
    reads         = ch_reads_filtered
    versions      = ch_versions
    multiqc_files = ch_multiqc_files
}
