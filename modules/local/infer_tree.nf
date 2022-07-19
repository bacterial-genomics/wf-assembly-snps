process INFER_TREE {

    publishDir "${params.outpath}/${prefix}",
        mode: "${params.publish_dir_mode}",
        pattern: "*.tree"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}" }

    container "staphb/parsnp@sha256:4f9ced31c7b7a4ef25046e4904c82d5489414f4ee5ce97e0a676788ea656c6df"

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