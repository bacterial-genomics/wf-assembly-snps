process CONVERT_GINGR_TO_FASTA_HARVESTTOOLS {

    tag { "${meta.snp_package}" }
    label "process_low"
    container "quay.io/biocontainers/parsnp@sha256:b46999fb9842f183443dd6226b461c1d8074d4c1391c1f2b1e51cc20cee8f8b2"

    input:
    tuple val(meta), path(snp_files)

    output:
    tuple val(meta), path("${meta.snp_package}.Core_Alignment.fasta")   , emit: core_alignment
    tuple val(meta), path("${meta.snp_package}.Gingr_to_FastA_File.tsv"), emit: qc_filecheck
    path(".command.{out,err}")
    path("versions.yml")                                            , emit: versions

    shell:
    '''
    source bash_functions.sh

    msg "INFO: Generating FastA file from Parsnp Gingr output file."

    harvesttools \
      -i "!{meta.snp_package}.ggr" \
      -M "!{meta.snp_package}.Core_Alignment.fasta"

    # Remove the 1 additional suffix Parsnp adds to the reference sample `.ref`
    if [[ $(grep -o -n '.ref' !{meta.snp_package}.Core_Alignment.fasta | wc -l) -eq 1 ]]; then
      sed -i 's/.ref//1' !{meta.snp_package}.Core_Alignment.fasta
    fi

    # Verify output
    echo -e "Sample name\tQC step\tOutcome (Pass/Fail)" > "!{meta.snp_package}.Gingr_to_FastA_File.tsv"
    if verify_minimum_file_size \
      "!{meta.snp_package}.Core_Alignment.fasta" \
      "!{meta.snp_package} Gingr to FastA File" \
      "!{params.min_gingr_to_fasta_filesize}"; then

      echo -e "NaN\t!{meta.snp_package} Gingr to FastA File\tPASS" >> "!{meta.snp_package}.Gingr_to_FastA_File.tsv"

    else
      echo -e "NaN\t!{meta.snp_package} Gingr to FastA File\tFAIL" >> "!{meta.snp_package}.Gingr_to_FastA_File.tsv"
    fi

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        harvesttools: $(harvesttools --version)
    END_VERSIONS
    '''
}
