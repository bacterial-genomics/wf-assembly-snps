process CALCULATE_PAIRWISE_DISTANCES_SNP_DISTS {

    tag ( "${meta.aligner}" )
    label "process_medium"
    container "staphb/snp-dists@sha256:9d457181cccacbbea0a3cb935edb1d066011fbc2f578694b6c5f9d9d58dcac15"

    input:
    tuple val(meta), path(snp_alignment)

    output:
    tuple val(meta), path("${meta.aligner}.SNP-Distances.Pairs.tsv"), emit: snp_distances
    path(".command.{out,err}")
    path("versions.yml")                                            , emit: versions

    shell:
    '''
    source bash_functions.sh

    msg "INFO: Calculating pairwise SNP distances."

    snp-dists \
      -m \
      -j !{task.cpus} \
      "!{meta.aligner}.SNPs.fa.gz" \
      | sort -nk3 \
      > "!{meta.aligner}.SNP-Distances.Pairs.tsv"

    # Reorder output and place self-vs-self matches at bottom of output file
    awk -F '\t' '$1 == $2 {print $0}' "!{meta.aligner}.SNP-Distances.Pairs.tsv" > self-pairs.tsv
    awk -F '\t' '$1 != $2 {print $0}' "!{meta.aligner}.SNP-Distances.Pairs.tsv" > non-self-pairs.tsv

    cat non-self-pairs.tsv self-pairs.tsv > "!{meta.aligner}.SNP-Distances.Pairs.tsv"

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        snp-dists: $(snp-dists -v | awk '{print $2}')
    END_VERSIONS
    '''
}
