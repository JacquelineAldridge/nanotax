include { SAMTOOLS_VIEW as BASECALL_FILTER     } from '../../../modules/nf-core/samtools/view/main'
include { PIGZ_COMPRESS as COMPRESS_CLASSIFIED } from '../../../modules/nf-core/pigz/compress/main'

include { DORADO_BASECALLER                    } from '../../../modules/local/dorado/basecaller/main'
include { DORADO_DEMUX                         } from '../../../modules/local/dorado/demux/main'


workflow BASECALLING {
    take:
    ch_pod5_dir              // channel: [ val(meta), path(pod5_dir) ]
    ch_samples               // channel: [ val(meta), path(fastq) ]
    val_dorado_barcoding_kit // string: dorado barcoding kit name

    main:
    ch_versions = channel.empty()

    DORADO_BASECALLER(ch_pod5_dir)

    BASECALL_FILTER(
        DORADO_BASECALLER.out.reads.map { meta, reads -> [meta, reads, []] },
        [[], []],
        [],
        [],
    )

    val_sample_sheet = ch_samples
        .collectFile(name: 'sample_sheet.csv', keepHeader: true) { meta, _fastq ->
            [
                "alias,barcode,kit,experiment_id,position_id\n",
                "${meta.id},${meta.barcode},${val_dorado_barcoding_kit},,\n",
            ].join('')
        }
        .first()

    DORADO_DEMUX(
        DORADO_BASECALLER.out.reads,
        val_sample_sheet,
        val_dorado_barcoding_kit,
    )


    ch_samples_with_sequences = DORADO_DEMUX.out.classified
        .flatMap { _meta, fastqs ->
            fastqs.collect { fastq ->
                [fastq.getSimpleName().replaceFirst(/^[^_]+_/, ""), fastq]
            }
        }
        .join(ch_samples.map { meta, _fastq -> [meta.id, meta] })
        .map { _id, fastq, meta -> [meta, fastq] }

    COMPRESS_CLASSIFIED(ch_samples_with_sequences)


    ch_versions = ch_versions.mix(
        DORADO_BASECALLER.out.versions,
        BASECALL_FILTER.out.versions,
        DORADO_DEMUX.out.versions,
    )

    emit:
    samples  = COMPRESS_CLASSIFIED.out.archive
    versions = ch_versions
}
