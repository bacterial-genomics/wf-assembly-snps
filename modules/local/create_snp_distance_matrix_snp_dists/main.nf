process CREATE_SNP_DISTANCE_MATRIX_SNP_DISTS {

    tag { "${meta.snp_package}" }
    label "process_low"
    container "staphb/snp-dists@sha256:9d457181cccacbbea0a3cb935edb1d066011fbc2f578694b6c5f9d9d58dcac15"

    input:
    tuple val(meta), path(snp_alignment_files)

    output:
    tuple val(meta), path ("*.SNP_Distances_Matrix.tsv"), emit: distance_matrix
    path (".command.{out,err}")
    path ("versions.yml")                               , emit: versions

    shell:
    '''
    source bash_functions.sh

    msg "INFO: Creating SNP distance matrix."

    snp-dists \
      -b \
      -j !{task.cpus} \
      "!{meta.snp_package}.SNPs.fa.gz" \
      > "!{meta.snp_package}.SNP_Distances_Matrix.tsv"

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        snp-dists: $(snp-dists -v | awk '{print $2}')
    END_VERSIONS
    '''
}
