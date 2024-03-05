process RECOMBINATION_GUBBINS {

    container "snads/gubbins@sha256:391a980312096f96d976f4be668d4dea7dda13115db004a50e49762accc0ec62"

    input:
    tuple val(meta_alignment), path(alignment)
    tuple val(meta_tree), path(starter_tree)

    output:
    tuple val(meta_tree), path("*.{gff,tree}"), emit: positions_and_tree
    path(".command.{out,err}")
    path("versions.yml")                      , emit: versions

    shell:
    '''
    source bash_functions.sh

    msg "INFO: Performing recombination using Gubbins."
    run_gubbins.py --starting-tree !{starter_tree} --prefix Gubbins !{alignment}

    # Rename output files
    mv Gubbins.recombination_predictions.gff Gubbins.recombination_positions.gff
    mv Gubbins.node_labelled.final_tree.tre Gubbins.labelled_tree.tree

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        gubbins: $(run_gubbins.py --version | sed 's/^/    /')
    END_VERSIONS
    '''
}
