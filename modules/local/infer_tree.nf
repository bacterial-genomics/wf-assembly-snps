process INFER_TREE {

    publishDir "${params.outpath}/${prefix}",
        mode: "${params.publish_dir_mode}",
        pattern: "*.tree"
    publishDir "${params.logpath}/command_outputs",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}" }

    container "snads/parsnp@sha256:0dc552de1cf91fb30aa25eb09b4d9eef841abae989760c937b0006dacd165377"

    input:
        path alignment
        val  prefix
        val  tree_method

    output:
        path("${prefix}.tree"), emit: tree
        path(".command.out")
        path(".command.err")
        path "versions.yml",    emit: versions

    shell:
    '''
    if [ "!{tree_method}" = "fasttree" ]; then
        fasttree -nt !{alignment} > !{prefix}.tree
    elif [ "!{tree_method}" = "raxml" ]; then
        raxmlHPC-PTHREADS -s !{alignment} -m GTRGAMMA -w !{workDir} -n !{prefix} -p 5280
        mv !{workDir}/RAxML_bestTree.!{prefix} !{prefix}.tree
    fi
    if [ ! -s !{prefix}.tree ]; then
        echo "Empty tree file. Did tree building program run out of RAM?" >> !{params.logpath}/stderr.nextflow.txt
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        fasttree: $(fasttree -expert &> tmp.txt; head -1 tmp.txt | sed 's/^/    /')
        raxml: $(raxmlHPC-PTHREADS -v | sed 's/^/    /')
    END_VERSIONS
    '''

}