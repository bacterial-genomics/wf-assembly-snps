process MASK_RECOMBINANT_POSITIONS_BIOPYTHON {

    tag( "${meta.snp_package}-${meta.recombination}" )
    label "process_medium"
    container "quay.io/biocontainers/biopython@sha256:10d755c731c82a22d91fc346f338ba47d5fd4f3b357828f5bbc903c9be865614"

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
      --rec_positions "!{meta.snp_package}-!{meta.recombination}.recombination_positions.txt" \
      --tree "!{meta.snp_package}-!{meta.recombination}.labelled_tree.tree" \
      > "!{meta.snp_package}-!{meta.recombination}.masked_recombination.fasta"

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        python: $(python --version)
    END_VERSIONS
    '''
}
