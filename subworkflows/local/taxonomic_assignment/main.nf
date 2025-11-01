include { EMU_ABUNDANCE } from '../../../modules/nf-core/emu/abundance/main'

workflow TAXONOMIC_ASSIGNMENT {
    take:
    ch_reads       // channel: [ val(meta), path(fastq) ]
    emu_database
    mmseqs2_database

    main:
    /*
     * EMU
     */
    EMU_ABUNDANCE(ch_reads, emu_database)

}
