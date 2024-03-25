process RECOMBINATION_CLONALFRAMEML {

    tag { "${meta.snp_package}" }
    label "process_medium"
    container "snads/clonalframeml@sha256:bc00db247195fdc6151793712a74cc9b272dc2c9f153bb0037415e387f15351e"

    input:
    tuple val(meta), path(core_alignment_fasta)
    tuple val(meta_alignment), path(alignment_files)

    output:
    tuple val(meta), path("*_{positions,tree}.*"), emit: positions_and_tree
    path(".command.{out,err}")
    path("versions.yml")                         , emit: versions

    shell:
    '''
    # ClonalFrameML needs tree labels to not be surrounded by single quotes
    sed -i "s/'//g" "!{meta.snp_package}.tree"

    ClonalFrameML \
      "!{meta.snp_package}.tree" \
      "!{core_alignment_fasta}" \
      "!{meta.snp_package}-ClonalFrameML"

    # Rename output file
    mv \
      "!{meta.snp_package}-ClonalFrameML.importation_status.txt" \
      "!{meta.snp_package}-ClonalFrameML.recombination_positions.txt"

    mv \
      "!{meta.snp_package}-ClonalFrameML.labelled_tree.newick" \
      "!{meta.snp_package}-ClonalFrameML.labelled_tree.tree"

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        clonalframeml: $(ClonalFrameML -version | sed 's/^/    /')
    END_VERSIONS
    '''
}
