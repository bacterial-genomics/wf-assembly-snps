process CONVERT_GINGR_TO_FASTA_HARVESTTOOLS {

    tag ( "${meta.aligner}" )
    container "quay.io/biocontainers/parsnp@sha256:b46999fb9842f183443dd6226b461c1d8074d4c1391c1f2b1e51cc20cee8f8b2"

    input:
    tuple val(meta), path(snp_files)

    output:
    tuple val(meta), path("${meta.aligner}.Core_Alignment.fasta")   , emit: core_alignment
    tuple val(meta), path("${meta.aligner}.Gingr_to_FastA_File.tsv"), emit: qc_filecheck
    path(".command.{out,err}")
    path("versions.yml")                                            , emit: versions

    shell:
    '''
    source bash_functions.sh

    msg "INFO: Generating FastA file from Parsnp Gingr output file."

    harvesttools \
      -i "!{meta.aligner}.ggr" \
      -M "!{meta.aligner}.Core_Alignment.fasta"

    # Remove all copies of `.ref` from FastA file
    sed -i 's/.ref//g' !{meta.aligner}.Core_Alignment.fasta

    # Verify output
    echo -e "Sample name\tQC step\tOutcome (Pass/Fail)" > "!{meta.aligner}.Gingr_to_FastA_File.tsv"
    if verify_minimum_file_size \
      "!{meta.aligner}.Core_Alignment.fasta" \
      "!{meta.aligner} Gingr to FastA File" \
      "!{params.min_gingr_to_fasta_filesize}"; then

      echo -e "NaN\t!{meta.aligner} Gingr to FastA File\tPASS" >> "!{meta.aligner}.Gingr_to_FastA_File.tsv"

    else
      echo -e "NaN\t!{meta.aligner} Gingr to FastA File\tFAIL" >> "!{meta.aligner}.Gingr_to_FastA_File.tsv"
    fi

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        harvesttools: $(harvesttools --version)
    END_VERSIONS
    '''
}
