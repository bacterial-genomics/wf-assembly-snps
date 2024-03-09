process MASK_RECOMBINANT_POSITIONS_BIOPYTHON {

    tag( "${meta.aligner}-${meta.recombination}" )
    container "gregorysprenger/biopython@sha256:77a50d5d901709923936af92a0b141d22867e3556ef4a99c7009a5e7e0101cc1"

    input:
    tuple val(meta), path(recombination_files)
    tuple val(meta_core_alignment), path(core_alignment_fasta)

    output:
    tuple val(meta), path("*.masked_recombination.fasta"), emit: masked_alignment
    path(".command.{out,err}")
    path("versions.yml")                                 , emit: versions

    shell:
    format = meta.recombination.toString().toLowerCase()
    '''
    source bash_functions.sh

    msg "INFO: Masking recombinant positions."

    mask_recombination.py \
      --alignment "!{core_alignment_fasta}" \
      --format !{format} \
      --rec_positions !{meta.recombination}.recombination_positions.* \
      --tree !{meta.recombination}.labelled_tree.* \
      > "!{meta.aligner}-!{meta.recombination}.masked_recombination.fasta"

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        python: $(python --version)
    END_VERSIONS
    '''
}
