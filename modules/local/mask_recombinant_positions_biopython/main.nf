process MASK_RECOMBINANT_POSITIONS_BIOPYTHON {

    tag( "${meta.aligner}-${meta.recombination}" )
    container "snads/mask-recombination@sha256:0df4f5e26b2beeb2a180c2e4d75148cde55d4cc62585b5053d6606c6623d33e4"

    input:
    tuple val(meta), path(recombination_files)
    tuple val(meta_alignment), path(alignment)

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
      --alignment !{alignment} \
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
