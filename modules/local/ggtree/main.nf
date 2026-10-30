process GGTREE {
    label 'process_single'

    container "https://depot.galaxyproject.org/singularity/bioconductor-ggtree:3.14.0--r44hdfd78af_0"

    input:
    tuple path(tree), path(clusters)

    output:
    path("tree.png"), emit: tree

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    tree.R ${tree} ${clusters}
    """
}
