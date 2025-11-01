include { OSFCLIENT_FETCH } from '../../../modules/nf-core/osfclient/fetch/main'

workflow PREPARE_DATABASES {
    take:
        val_skip_emu
        emu_db_name

    main:
        ch_versions = channel.empty()
        ch_emu_db = channel.empty()

        if (!val_skip_emu) {
            if (emu_db_name == 'default') {
                OSFCLIENT_FETCH([
                    [id: 'emu_database'],
                    'g6w5e',
                    '/database/emu/emu_default_db.tar.gz'
                ])
                ch_versions = ch_versions.mix(OSFCLIENT_FETCH.out.versions)

                ch_emu_db = OSFCLIENT_FETCH.out.download_files
            } else {
                error "Currently only 'default' emu database is supported"
            }
        }


    emit:
        emu_database = ch_emu_db
}
