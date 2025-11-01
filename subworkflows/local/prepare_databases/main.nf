include { OSFCLIENT_FETCH } from '../../../modules/nf-core/osfclient/fetch/main'
include { UNTAR           } from '../../../modules/nf-core/untar/main'

workflow PREPARE_DATABASES {
    take:
    val_skip_emu
    emu_db_name

    main:
    ch_versions = channel.empty()
    ch_emu_db = channel.empty()

    if (!val_skip_emu) {
        db_paths = [
            'default': 'osfstorage/emu-prebuilt/emu.tar.gz',
            'rdp': 'osfstorage/emu-prebuilt/rdp.tar.gz',
            'silva': 'osfstorage/emu-prebuilt/silva_database.tar.gz',
        ]

        OSFCLIENT_FETCH(
            [
                [id: "emu_database_${emu_db_name}"],
                '56uf7',
                db_paths[emu_db_name],
            ]
        )
        ch_versions = ch_versions.mix(OSFCLIENT_FETCH.out.versions)

        UNTAR(OSFCLIENT_FETCH.out.download_files)

        ch_emu_db = UNTAR.out.untar.map { _meta, file -> file }
    }

    emit:
    emu_database = ch_emu_db
}
