include { EMU_ABUNDANCE     } from '../../../modules/nf-core/emu/abundance/main'
include { MMSEQS_EASYSEARCH } from '../../../modules/nf-core/mmseqs/easysearch/main'

workflow TAXONOMIC_ASSIGNMENT {
    take:
    ch_reads // channel: [ val(meta), path(fastq) ]
    ch_emu_db
    ch_mmseqs2_db
    val_skip_emu
    val_skip_mmseqs2

    main:
    ch_versions = channel.empty()

    /*
     * EMU
     */
    if (!val_skip_emu) {
        EMU_ABUNDANCE(ch_reads, ch_emu_db.map { _meta, file -> file })
        ch_versions = ch_versions.mix(EMU_ABUNDANCE.out.versions)
    }

    if (!val_skip_mmseqs2) {
        MMSEQS_EASYSEARCH(ch_reads, ch_mmseqs2_db)
        ch_versions = ch_versions.mix(MMSEQS_EASYSEARCH.out.versions)
    }

    emit:
    versions = ch_versions
}
