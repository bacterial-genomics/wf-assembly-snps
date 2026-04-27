process CLUSTER {
    label 'process_single'

    container "https://depot.galaxyproject.org/singularity/numpy%3A2.2.2"

    input:
    path(snpmatrix)
    val(threshold)

    output:
    path("clusters.txt"), emit: clusters

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    outbreak_detection.py -i ${snpmatrix} -d ${threshold} > clusters.txt
    """
}
