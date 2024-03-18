process CREATE_MASKED_SNP_DISTANCE_MATRIX_SNP_DISTS {

    tag { "${meta.snp_package}-${meta.recombination}" }
    label "process_medium"
    container "staphb/snp-dists@sha256:9d457181cccacbbea0a3cb935edb1d066011fbc2f578694b6c5f9d9d58dcac15"

    input:
    tuple val(meta), path(masked_alignment)

    output:
    tuple val(meta), path ("*.Masked_SNP_Distances_Matrix.tsv"), emit: distance_matrix
    path (".command.{out,err}")
    path ("versions.yml")                                      , emit: versions

    shell:
    '''
    source bash_functions.sh

    msg "INFO: Creating SNP distance matrix using the masked multi-fasta file."

    snp-dists \
      -b \
      -k \
      -j !{task.cpus} \
      "!{masked_alignment}" \
      > "!{meta.snp_package}-!{meta.recombination}.Masked_SNP_Distances_Matrix.tsv"

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        snp-dists: $(snp-dists -v | awk '{print $2}')
    END_VERSIONS
    '''
}
