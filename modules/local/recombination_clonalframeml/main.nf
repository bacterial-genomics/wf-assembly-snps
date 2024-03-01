process RECOMBINATION_CLONALFRAMEML {

    container "snads/clonalframeml@sha256:bc00db247195fdc6151793712a74cc9b272dc2c9f153bb0037415e387f15351e"

    input:
    tuple val(meta_alignment), path(alignment)
    tuple val(meta_tree), path(starter_tree)

    output:
    tuple val(meta_tree), path("*_{positions,tree}.*"), emit: positions_and_tree
    path(".command.{out,err}")
    path("versions.yml")                              , emit: versions

    shell:
    '''
    # ClonalFrameML needs tree labels to not be surrounded by single quotes
    sed -i "s/'//g" !{starter_tree}

    ClonalFrameML !{starter_tree} !{alignment} ClonalFrameML

    # Rename output file
    mv ClonalFrameML.importation_status.txt ClonalFrameML.recombination_positions.txt

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        clonalframeml: $(ClonalFrameML -version | sed 's/^/    /')
    END_VERSIONS
    '''
}
