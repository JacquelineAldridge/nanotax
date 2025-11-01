process MMSEQS_CREATETAXDB {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mmseqs2:17.b804f--hd6d6fdc_1'
        : 'biocontainers/mmseqs2:17.b804f--hd6d6fdc_1'}"

    input:
    tuple val(meta), path(db)
    tuple val(meta2), path(taxdump_dir)
    tuple val(meta3), path(tax_mapping_file)

    output:
    tuple val(meta), path(db), emit: db_with_taxonomy
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: "*.dbtype"
    taxdump_opt = taxdump_dir ? "--ncbi-tax-dump ${taxdump_dir}" : ""
    tax_mapping_opt = taxdump_dir && tax_mapping_file ? "--tax-mapping-file ${tax_mapping_file}" : ""
    """
    DB_INPUT_PATH_NAME=\$(find -L "${db}/" -maxdepth 1 -name "${args2}" | sed 's/\\.[^.]*\$//' |  sed -e 'N;s/^\\(.*\\).*\\n\\1.*\$/\\1\\n\\1/;D' )

    mmseqs createtaxdb \\
      \${DB_INPUT_PATH_NAME} \\
      ./tmp \\
      --threads ${task.cpus} \\
      ${taxdump_opt} \\
      ${tax_mapping_opt} \\
      ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """
}
