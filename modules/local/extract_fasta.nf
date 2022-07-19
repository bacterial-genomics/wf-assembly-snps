process EXTRACT_FASTA {

    publishDir "${params.outpath}/parsnp",
        mode: "${params.publish_dir_mode}",
        pattern: "*.fasta"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}" }

    container "snads/xmfa-to-fasta@sha256:ef61b6d2c1a3ac675ecd102d64488f65f715745788deaf9fae2d7ab69c71c277"

    input:
        path xmfa

    output:
        path "parsnp.fasta", emit: parsnp_fasta
        path ".command.out"
        path ".command.err"
        path "versions.yml", emit: versions

    shell:
    '''
    convert_xmfa_to_fasta.py --xmfa !{xmfa} > parsnp.fasta
    sed -i.bak 's/.ref//g' parsnp.fasta
    rm parsnp.fasta.bak

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        python: $(python --version)
    END_VERSIONS
    '''

}