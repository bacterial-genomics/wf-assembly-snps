process BUILD_PHYLOGENETIC_TREE_PARSNP {

    container "staphb/parsnp@sha256:4f9ced31c7b7a4ef25046e4904c82d5489414f4ee5ce97e0a676788ea656c6df"

    input:
    tuple val(method), path(masked_alignment)

    output:
    tuple val(method), path("${method}.Tree_Output_File.tsv"), emit: qc_filecheck
    tuple val(method), path("${method}.tree")                , emit: tree
    path(".command.{out,err}")
    path("versions.yml")                                     , emit: versions

    shell:
    '''
    source bash_functions.sh

    if [[ "!{params.tree_method}" = "fasttree" ]]; then
      msg "INFO: Building phylogenetic tree using FastTree."

      fasttree \
        -nt !{masked_alignment} \
        > !{method}.tree

    elif [[ "!{params.tree_method}" = "raxml" ]]; then
      msg "INFO: Building phylogenetic tree using RaxML."

      raxmlHPC-PTHREADS \
        -s !{masked_alignment} \
        -m GTRGAMMA \
        -n !{method} \
        -p 5280

      mv RAxML_bestTree.!{method} !{method}.tree
    fi

    echo -e "Sample name\tQC step\tOutcome (Pass/Fail)" > "!{method}.Tree_Output_File.tsv"
    if verify_minimum_file_size "!{method}.tree" "Final !{method} Tree Output" "!{params.min_tree_filesize}"; then
      echo -e "NaN\tFinal !{method} Tree Output\tPASS" >> "!{method}.Tree_Output_File.tsv"
    else
      echo -e "NaN\tFinal !{method} Tree Output\tFAIL" >> "!{method}.Tree_Output_File.tsv"
    fi

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        fasttree: $(fasttree -expert &> tmp.txt; head -1 tmp.txt | sed 's/^/    /')
        raxml: $(raxmlHPC-PTHREADS -v | sed 's/^/    /')
    END_VERSIONS
    '''

}