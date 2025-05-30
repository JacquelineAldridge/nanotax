include { DORADO_BASECALLER                } from '../../modules/local/dorado/basecaller/main.nf'
include { DORADO_DEMUX                     } from '../../modules/local/dorado/demux/main.nf'
include { SAMTOOLS_VIEW as BASECALL_FILTER } from '../../modules/nf-core/samtools/view/main.nf'

workflow BASECALLING {
    take:
    ch_pod5_dir // pod5 directory from sequencing
    ch_samples  // channel: [ val(meta)]

    main:
    ch_versions = Channel.empty()

    DORADO_BASECALLER(ch_pod5_dir)
    BASECALL_FILTER(
        DORADO_BASECALLER.out.reads.map { meta, reads -> [meta, reads, null] },
        [[:], null],
        [],
        'csi',
    )


    ch_versions = ch_versions.mix(
        DORADO_BASECALLER.out.versions,
        BASECALL_FILTER.out.versions,
    )

    emit:
    versions = ch_versions
}
