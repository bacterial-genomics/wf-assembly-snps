process CALCULATE_PAIRWISE_DISTANCES_BIOPYTHON {

    tag ( "${meta.aligner}" )
    label "process_medium"
    container "snads/hamming-dist@sha256:3ecbf4f963adfd8de843f57487ec68ed71614d62956ce4993af3679d08785c48"

    input:
    tuple val(meta), path(snp_files)

    output:
    tuple val(meta), path("${meta.aligner}.Pairwise_SNP_Distances_File.tsv"), emit: qc_filecheck
    tuple val(meta), path("${meta.aligner}.SNP-Distances.Pairs.tsv")        , emit: snp_distances
    path(".command.{out,err}")
    path("versions.yml")                                                    , emit: versions

    shell:
    '''
    source bash_functions.sh

    msg "INFO: Calculating pairwise SNP distances from all sample pairs."

    pairwiseDistances.py \
      -n "!{task.cpus}" \
      "!{meta.aligner}.SNPs.fa" \
      | sort -k3,3n \
      > "!{meta.aligner}.SNP-Distances.Pairs.tsv"

    echo -e "Sample name\tQC step\tOutcome (Pass/Fail)" > "!{meta.aligner}.Pairwise_SNP_Distances_File.tsv"
    if verify_minimum_file_size "!{meta.aligner}.SNP-Distances.Pairs.tsv" 'Pairwise SNP Distances' "!{params.min_snp_distance_filesize}"; then
      echo -e "!{meta.aligner}\tPairwise SNP Distances\tPASS" >> "!{meta.aligner}.Pairwise_SNP_Distances_File.tsv"
    else
      echo -e "!{meta.aligner}\tPairwise SNP Distances\tFAIL" >> "!{meta.aligner}.Pairwise_SNP_Distances_File.tsv"
    fi

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        python: $(python --version)
    END_VERSIONS
    '''
}
