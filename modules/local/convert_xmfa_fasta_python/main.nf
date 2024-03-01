process CONVERT_XMFA_FASTA_PYTHON {

    container "snads/xmfa-to-fasta@sha256:ef61b6d2c1a3ac675ecd102d64488f65f715745788deaf9fae2d7ab69c71c277"

    input:
    tuple val(meta), path(xmfa)

    output:
    tuple val(meta), path("${meta.aligner}.Core_Alignment.fasta"), emit: core_alignment
    path(".command.{out,err}")
    path("versions.yml")                                         , emit: versions

    shell:
    '''
    source bash_functions.sh

    msg "INFO: Generating FastA file from Parsnp XMFA output file."

    convert_xmfa_to_fasta.py \
      --xmfa !{xmfa} \
      > !{meta.aligner}.Core_Alignment.fasta

    # Remove all copies of `.ref` from FastA file
    sed -i 's/.ref//g' !{meta.aligner}.Core_Alignment.fasta

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        python: $(python --version)
    END_VERSIONS
    '''
}
