process BUILD_PHYLOGENETIC_TREE_PARSNP {

    container "quay.io/biocontainers/parsnp@sha256:b46999fb9842f183443dd6226b461c1d8074d4c1391c1f2b1e51cc20cee8f8b2"

    input:
    tuple val(meta), path(masked_alignment)

    output:
    tuple val(meta), path("${meta.recombination}.Tree_Output_File.tsv"), emit: qc_filecheck
    tuple val(meta), path("${meta.recombination}.Final.tree")          , emit: tree
    path(".command.{out,err}")
    path("versions.yml")                                               , emit: versions

    shell:
    '''
    source bash_functions.sh

    if [[ "!{params.tree_method}" = "fasttree" ]]; then
      msg "INFO: Building phylogenetic tree using FastTree."

      fasttree \
        -nt !{masked_alignment} \
        > !{meta.recombination}.Final.tree

    elif [[ "!{params.tree_method}" = "raxml" ]]; then
      msg "INFO: Building phylogenetic tree using RaxML."

      raxmlHPC-PTHREADS \
        -s !{masked_alignment} \
        -m GTRGAMMA \
        -n !{meta.recombination} \
        -p 5280

      mv RAxML_bestTree.!{meta.recombination} !{meta.recombination}.Final.tree
    fi

    echo -e "Sample name\tQC step\tOutcome (Pass/Fail)" > "!{meta.recombination}.Tree_Output_File.tsv"
    if verify_minimum_file_size "!{meta.recombination}.tree" "Final !{meta.recombination} Tree Output" "!{params.min_tree_filesize}"; then
      echo -e "!{meta.recombination}\tFinal Tree Output\tPASS" >> "!{meta.recombination}.Tree_Output_File.tsv"
    else
      echo -e "!{meta.recombination}\tFinal Tree Output\tFAIL" >> "!{meta.recombination}.Tree_Output_File.tsv"
    fi

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        fasttree: $(fasttree -expert &> tmp.txt; head -1 tmp.txt | sed 's/^/    /')
        raxml: $(raxmlHPC-PTHREADS -v | sed 's/^/    /')
    END_VERSIONS
    '''
}
