process EXTRACT_SNP_POSITIONS_PARSNP {

    container "staphb/parsnp@sha256:4f9ced31c7b7a4ef25046e4904c82d5489414f4ee5ce97e0a676788ea656c6df"

    input:
    tuple val(meta), path(gingr_alignment)

    output:
    tuple val(meta), path("Extracted_SNP_Positions.tsv"), emit: qc_filecheck
    tuple val(meta), path("Parsnp.SNPs.fa")             , emit: snps
    path(".command.{out,err}")
    path("versions.yml")                                , emit: versions

    shell:
    '''
    source bash_functions.sh

    msg "INFO: Extracting SNP positions."

    harvesttools \
      -i "!{gingr_alignment}" \
      -S Parsnp.SNPs.fa

    echo -e "Sample name\tQC step\tOutcome (Pass/Fail)" > Extracted_SNP_Positions.tsv
    if verify_minimum_file_size "Parsnp.SNPs.fa" 'Extracted SNP Positions' "!{params.min_snp_positions_filesize}"; then
      echo -e "NaN\tExtracted SNP Positions\tPASS" >> Extracted_SNP_Positions.tsv
    else
      echo -e "NaN\tExtracted SNP Positions\tFAIL" >> Extracted_SNP_Positions.tsv
    fi

    # Remove '.ref' suffix from reference genome in SNPs file
    sed -i "s/\\.ref//1" Parsnp.SNPs.fa

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        harvesttools: $(harvesttools --version)
    END_VERSIONS
    '''

}
