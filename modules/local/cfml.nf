process CFML {

    publishDir "${params.outpath}/clonalframeml",
        mode: "${params.publish_dir_mode}",
        pattern: "clonalframeml.*"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}" }

    container "snads/clonalframeml@sha256:bc00db247195fdc6151793712a74cc9b272dc2c9f153bb0037415e387f15351e"

    input:
        path alignment
        path starter_tree

    output:
        path "clonalframeml.importation_status.txt", emit: recombination_positions
        path "clonalframeml.labelled_tree.newick",   emit: node_labelled_tree
        path ".command.out"
        path ".command.err"
        path "versions.yml",                         emit: versions

    shell:
    '''
    # ClonalFrameML needs tree labels to not be surrounded by single quotes
    sed -i.bak "s/'//g" !{starter_tree}
    rm !{starter_tree}.bak
    ClonalFrameML !{starter_tree} !{alignment} clonalframeml

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        clonalframeml: $(ClonalFrameML -version | sed 's/^/    /')
    END_VERSIONS
    '''

}