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
    container "staphb/parsnp@sha256:4f9ced31c7b7a4ef25046e4904c82d5489414f4ee5ce97e0a676788ea656c6df"

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
