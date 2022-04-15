process DISTANCE_MATRIX {

    publishDir "${params.outpath}/parsnp",
        mode: "${params.publish_dir_mode}",
        pattern: "*.tsv"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}" }

    params.enable_conda_yml ? "$baseDir/conda/linux/python3.yml" : null
    // conda 'conda-forge::python=3.10.1'
    container "snads/hamming-dist@sha256:3ecbf4f963adfd8de843f57487ec68ed71614d62956ce4993af3679d08785c48"

    input:
        path snp_distances

    output:
        path "SNP-distances.matrix.tsv"
        path ".command.out"
        path ".command.err"
        path "versions.yml", emit: versions

    shell:
    '''
    source bash_functions.sh

    # Create a 2-dimensional matrix from pairwise SNP distances of all sample pairs
    snp_matrix=SNP-distances.matrix.tsv

    msg "INFO: starting to form SNP matrix"
    pairwiseTo2d.py \
     -i "!{snp_distances}" \
     -o "${snp_matrix}" \
     --sort
    msg "INFO: finished SNP matrix formation"

    if ! verify_file_minimum_size "${snp_matrix}" 'SNP matrix' '20c'; then
      msg "ERROR: SNP distance matrix not formed" >&2
      exit 1
    fi
    sed -i "s/\t-/\t0/g" "${snp_matrix}"

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        python: $(python --version)
    END_VERSIONS
    '''

}