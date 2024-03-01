process CALCULATE_PAIRWISE_DISTANCES_BIOPYTHON {

    label "process_medium"

    container "snads/hamming-dist@sha256:3ecbf4f963adfd8de843f57487ec68ed71614d62956ce4993af3679d08785c48"

    input:
    tuple val(meta), path(snps)

    output:
    tuple val(meta), path("Pairwise_SNP_Distances_File.tsv"), emit: qc_filecheck
    tuple val(meta), path("SNP-distances.pairs.tsv")        , emit: snp_distances
    path(".command.{out,err}")
    path("versions.yml")                                    , emit: versions

    shell:
    '''
    source bash_functions.sh

    msg "INFO: Calculating pairwise SNP distances from all sample pairs."

    pairwiseDistances.py \
      -n "!{task.cpus}" \
      "!{snps}" \
      | sort -k3,3n \
      > SNP-distances.pairs.tsv

    echo -e "Sample name\tQC step\tOutcome (Pass/Fail)" > Pairwise_SNP_Distances_File.tsv
    if verify_minimum_file_size "SNP-distances.pairs.tsv" 'Pairwise SNP Distances' "!{params.min_snp_distance_filesize}"; then
      echo -e "!{meta.aligner}\tPairwise SNP Distances\tPASS" >> Pairwise_SNP_Distances_File.tsv
    else
      echo -e "!{meta.aligner}\tPairwise SNP Distances\tFAIL" >> Pairwise_SNP_Distances_File.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        python: $(python --version)
    END_VERSIONS
    '''
}
