process GUBBINS {

    publishDir "${params.outpath}/gubbins",
        mode: "${params.publish_dir_mode}",
        pattern: "gubbins.*"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}" }

    container = "snads/gubbins@sha256:391a980312096f96d976f4be668d4dea7dda13115db004a50e49762accc0ec62"

    input:
        path alignment
        path starter_tree

    output:
        path "gubbins.recombination_predictions.gff", emit: recombination_positions
        path "gubbins.node_labelled.final_tree.tre",  emit: node_labelled_tree
        path ".command.out"
        path ".command.err"
        path "versions.yml",                          emit: versions

    shell:
    '''
    run_gubbins.py --starting-tree !{starter_tree} --prefix gubbins !{alignment}

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        gubbins: $(run_gubbins.py --version | sed 's/^/    /')
    END_VERSIONS
    '''

}