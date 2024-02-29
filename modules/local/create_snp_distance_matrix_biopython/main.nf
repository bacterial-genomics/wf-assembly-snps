process CREATE_SNP_DISTANCE_MATRIX_BIOPYTHON {

    container "snads/hamming-dist@sha256:3ecbf4f963adfd8de843f57487ec68ed71614d62956ce4993af3679d08785c48"

    input:
    path(snp_distances)

    output:
    path ("SNP-distances.matrix.tsv")
    path("${snp_distances}.gz")
    path (".command.{out,err}")
    path ("versions.yml")            , emit: versions

    shell:
    '''
    source bash_functions.sh

    msg "INFO: Using SNP distances to create a distance matrix."

    pairwiseTo2d.py \
      -i "!{snp_distances}" \
      -o SNP-distances.matrix.tsv \
      --sort

    echo -e "Sample name\tQC step\tOutcome (Pass/Fail)" > SNP_Distance_Matrix_File.tsv
    if verify_minimum_file_size "SNP-distances.matrix.tsv" 'SNP Distance Matrix' "!{params.min_distance_matrix_filesize}"; then
      echo -e "NaN\tSNP Distance Matrix\tPASS" >> SNP_Distance_Matrix_File.tsv
    else
      echo -e "NaN\tSNP Distance Matrix\tFAIL" >> SNP_Distance_Matrix_File.tsv
    fi

    sed -i "s/\t-/\t0/g" SNP-distances.matrix.tsv

    # Gzip compress SNP distances
    gzip -9f !{snp_distances}

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        python: $(python --version)
    END_VERSIONS
    '''
}
