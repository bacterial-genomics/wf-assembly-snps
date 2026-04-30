process PATRISTICDISTANCE {
    label 'process_low'

    container "https://depot.galaxyproject.org/singularity/biopython%3A1.84"

    input:
    path(newick)

    output:
    path("*.tsv")                , emit: tsv
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    patristic_distance.py \\
        -i $newick

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$(python --version 2>&1 | sed 's/Python //;')
    END_VERSIONS

    """
}
