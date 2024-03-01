process CORE_GENOME_ALIGNMENT_PARSNP {

    label "process_medium"
    container "staphb/parsnp@sha256:4f9ced31c7b7a4ef25046e4904c82d5489414f4ee5ce97e0a676788ea656c6df"

    input:
    tuple val(meta_input)    , path(input_files)
    tuple val(meta_reference), path(reference_file)

    output:
    tuple val(meta_input), path("Parsnp_Alignment_File.tsv"), emit: qc_filecheck
    tuple val(meta_input), path("Parsnp.ggr")               , emit: gingr_alignment
    tuple val(meta_input), path("Parsnp.xmfa")              , emit: core_alignment
    tuple val(meta_input), path("Parsnp.tree")              , emit: phylogeny
    path(".command.{out,err}")
    path("versions.yml")                                    , emit: versions

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
    echo -e "Sample name\tQC step\tOutcome (Pass/Fail)" > Parsnp_Alignment_File.tsv
    for file in Parsnp.ggr Parsnp.xmfa; do
      if [ ${file#*.} == "ggr" ]; then output="Gingr"; else output="Core"; fi

      if verify_minimum_file_size "${file}" "Parsnp ${output} Alignment File" "!{params.min_parsnp_alignment_filesize}"; then
        echo -e "NaN\t${output} Alignment File\tPASS" >> Parsnp_Alignment_File.tsv
      else
        echo -e "NaN\t${output} Alignment File\tFAIL" >> Parsnp_Alignment_File.tsv
      fi
    done

    # Remove reference label from tip in tree file
    sed -i "s/.ref//1" Parsnp.tree

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        parsnp: $(parsnp --version | sed 's/^/    /')
    END_VERSIONS
    '''
}
