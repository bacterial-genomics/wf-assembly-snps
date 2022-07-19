process PAIRWISE_DISTANCES {

    publishDir "${params.outpath}/parsnp",
        mode: "${params.publish_dir_mode}",
        pattern: "*.tsv"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}" }

    label "process_medium"
    label "error_retry"

    params.enable_conda_yml ? "$baseDir/conda/linux/NEEDS-NEWFILE.yml" : null
    // conda 'bioconda::FIXME'
    container "snads/hamming-dist@sha256:3ecbf4f963adfd8de843f57487ec68ed71614d62956ce4993af3679d08785c48"

    input:
        path snps_file

    output:
        path "SNP-distances.pairs.tsv", emit: snp_distances
        path ".command.out"
        path ".command.err"
        path "versions.yml", emit: versions

    shell:
    '''
    source bash_functions.sh

    # Calculate pairwise SNP distances from all sample pairs
    snp_distances="SNP-distances.pairs.tsv"
    msg "INFO: starting to calculate pairwise SNP distances"
    pairwiseDistances.py \
      -n "!{task.cpus}" \
      "!{snps_file}" \
       | sort -k3,3n \
       > "${snp_distances}"
    msg "INFO: finished calculating pairwise SNP distances"

    if ! verify_file_minimum_size "${snp_distances}" 'SNP distances' '20c'; then
      msg "ERROR: pairwise SNP distances not calculated" >&2
      exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        python: $(python --version)
    END_VERSIONS
    '''

}