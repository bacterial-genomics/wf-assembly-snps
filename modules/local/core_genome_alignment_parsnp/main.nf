process CORE_GENOME_ALIGNMENT_PARSNP {

    label "process_medium"

    container "staphb/parsnp@sha256:4f9ced31c7b7a4ef25046e4904c82d5489414f4ee5ce97e0a676788ea656c6df"

    input:
    path(input_files)
    path(reference_file)

    output:
    path("Parsnp_{Core,Gingr}_Alignment_File.tsv"), emit: qc_filecheck
    path("Parsnp.ggr")                            , emit: gingr_alignment
    path("Parsnp.xmfa")                           , emit: core_alignment
    path("Parsnp.tree")                           , emit: phylogeny
    path(".command.{out,err}")
    path("versions.yml")                          , emit: versions

    shell:
    curatedInput = params.curated_input             ? "-c"             : ""
    treeMethod   = params.tree_method == "fasttree" ? "--use-fasttree" : ""
    '''
    source bash_functions.sh

    # TEMPORARY hack to avoid the v1.2 issue that thinks input ref file has
    #  aligned nucleotide sequences if the deflines contain hyphens
    sed -i "s/-//g" !{reference_file}

    msg "INFO: Performing the parsnp system call."

    parsnp \
      -v \
      -d !{input_files} \
      -r !{reference_file} \
      -o Parsnp \
      -p !{task.cpus} \
      -P !{params.max_partition_size} \
      --verbose \
      !{treeMethod} \
      !{curatedInput}

    # Move files to current directory
    mv Parsnp/parsnp.ggr Parsnp.ggr
    mv Parsnp/parsnp.xmfa Parsnp.xmfa
    mv Parsnp/parsnp.tree Parsnp.tree

    # Verify output
    echo -e "Sample name\tQC step\tOutcome (Pass/Fail)" > Parsnp_Gingr_Alignment_File.tsv
    if verify_minimum_file_size "Parsnp.ggr" 'Parsnp Gingr Alignment File' "!{params.min_tree_filesize}"; then
      echo -e "NaN\tParsnp Gingr Alignment File\tPASS" >> Parsnp_Gingr_Alignment_File.tsv
    else
      echo -e "NaN\tParsnp Gingr Alignment File\tFAIL" >> Parsnp_Gingr_Alignment_File.tsv
    fi

    echo -e "Sample name\tQC step\tOutcome (Pass/Fail)" > Parsnp_Core_Alignment_File.tsv
    if verify_minimum_file_size "Parsnp.xmfa" 'Parsnp Core Alignment File' "!{params.min_core_alignment_filesize}"; then
      echo -e "NaN\tParsnp Gingr Alignment File\tPASS" >> Parsnp_Core_Alignment_File.tsv
    else
      echo -e "NaN\tParsnp Gingr Alignment File\tFAIL" >> Parsnp_Core_Alignment_File.tsv
    fi

    # Remove reference label from tip in tree file
    sed -i "s/.ref//1" Parsnp.tree

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        parsnp: $(parsnp --version | sed 's/^/    /')
    END_VERSIONS
    '''
}
