process EXTRACT_SNPS {

    publishDir "${params.outpath}/parsnp",
        mode: "${params.publish_dir_mode}",
        pattern: "*.fa"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}" }

    params.enable_conda_yml ? "$baseDir/conda/linux/harvesttools.yml" : null
    // conda 'bioconda::harvesttools=1.2'
    container "snads/parsnp@sha256:0dc552de1cf91fb30aa25eb09b4d9eef841abae989760c937b0006dacd165377"

    input:
        path ggr

    output:
        path "SNPs.fa",      emit: snps_file
        path ".command.out"
        path ".command.err"
        path "versions.yml", emit: versions

    shell:
    '''
    source bash_functions.sh

    # Extract only SNP positions from the gingr file
    msg "INFO: starting SNP extraction"
    harvesttools \
      -i "!{ggr}" \
      -S SNPs.fa
    msg "INFO: finished SNP extraction"

    if ! verify_file_minimum_size SNPs.fa "SNP" "10c"; then
      msg "ERROR: harvesttools failed to extract SNPs" >&2
      exit 1
    fi

    # Remove '.ref' suffix from reference genome in SNPs file
    sed -i "s/\\.ref//1" SNPs.fa

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        harvesttools: $(harvesttools --version)
    END_VERSIONS
    '''

}
