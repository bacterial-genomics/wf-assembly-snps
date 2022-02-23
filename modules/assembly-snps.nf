nextflow.enable.dsl = 2


process FIND_INFILES {

    input:
        path(inpath)

    output:
        path("find_infiles.success.txt"), emit: find_infiles_success

    script:
    """
    find_infiles.sh ${inpath}
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """

}


process INFILE_HANDLING {
    
    input:
        path(find_infiles_success)
        path(inpath)
        val refpath

    output:
        path(".ref"), emit: refpath
        path(".tmp"), emit: tmppath
        path(".tmp/*")

    script:
    """
    infile_handling.sh ${inpath} ${refpath}
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """
    
    stub:
    """
    if [[ ! -d .tmp || .ref ]]; then
        cp -r ${params.examplepath}/.tmp .
        cp -r ${params.examplepath}/.ref .
    fi
    """

}


process RUN_PARSNP {
    publishDir "${params.outpath}", mode: "copy"

    params.enable_conda_yml ? "$baseDir/conda/linux/parsnp.yml" : null
    //conda 'bioconda::parsnp=1.1.3'
    container "snads/parsnp@sha256:0dc552de1cf91fb30aa25eb09b4d9eef841abae989760c937b0006dacd165377"

    input:
        path(refpath)
        path(tmppath)

    output:
        path("parsnp/parsnp.ggr"), emit: parsnp_ggr
        path("parsnp/parsnp.xmfa"), emit: parsnp_xmfa
        path("parsnp/parsnp.tree"), emit: parsnp_tree

    script:
    """
    # Log version info
    echo "RUN_PARSNP:" >> ${params.logpath}/versions.txt
    parsnp --version | sed 's/^/    /' >> ${params.logpath}/versions.txt

    run_parsnp.sh ${tmppath} parsnp ${refpath}/* ${task.cpus}
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """

    stub:
    """
    if [[ ! -f parsnp/parsnp.ggr || parsnp/parsnp.tree || parsnp/parsnp.xmfa ]]; then
        mkdir parsnp
        cp ${params.examplepath}/parsnp/parsnp.ggr parsnp/
        cp ${params.examplepath}/parsnp/parsnp.xmfa parsnp/
        cp ${params.examplepath}/parsnp/parsnp.tree parsnp/
    fi
    """

}


process EXTRACT_SNPS {
    publishDir "${params.outpath}/parsnp", mode: "copy"

    params.enable_conda_yml ? "$baseDir/conda/linux/harvesttools.yml" : null
    // conda 'bioconda::harvesttools=1.2'
    container "snads/parsnp@sha256:0dc552de1cf91fb30aa25eb09b4d9eef841abae989760c937b0006dacd165377"

    input:
        path(parsnp_ggr)

    output:
        path("SNPs.fa"), emit: snps_file

    script:
    """
    # Log version info
    echo "EXTRACT_SNPS:" >> ${params.logpath}/versions.txt
    parsnp --version | sed 's/^/    /' >> ${params.logpath}/versions.txt

    extract_snps.sh "${parsnp_ggr}" "SNPs.fa"
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """

    stub:
    """
    if [[ ! -f SNPs.fa ]]; then
        cp ${params.examplepath}/parsnp/SNPs.fa .
    fi
    """

}


process PAIRWISE_DISTANCES {
    publishDir "${params.outpath}/parsnp", mode: "copy"

    params.enable_conda_yml ? "$baseDir/conda/linux/NEEDS-NEWFILE.yml" : null
    // conda 'bioconda::FIXME'
    container "snads/hamming-dist@sha256:3ecbf4f963adfd8de843f57487ec68ed71614d62956ce4993af3679d08785c48"

    input:
        path(snps_file)

    output:
        path("SNP-distances.pairs.tsv"), emit: snp_distances

    script:
    """
    pairwise_distances.sh ${task.cpus} ${snps_file}
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """

    stub:
    """
    if [[ ! -f SNP-distances.pairs.tsv ]]; then
        cp ${params.examplepath}/parsnp/SNP-distances.pairs.tsv .
    fi
    """

}


process DISTANCE_MATRIX {
    publishDir "${params.outpath}/parsnp", mode: "copy"

    params.enable_conda_yml ? "$baseDir/conda/linux/python3.yml" : null
    // conda 'conda-forge::python=3.10.1'
    container "snads/hamming-dist@sha256:3ecbf4f963adfd8de843f57487ec68ed71614d62956ce4993af3679d08785c48"

    input:
        path(snp_distances)

    output:
        path("SNP-distances.matrix.tsv")

    script:
    """
    distance_matrix.sh ${snp_distances}
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """

    stub:
    """
    if [[ ! -f SNP-distances.matrix.tsv ]]; then
        cp ${params.examplepath}/parsnp/SNP-distances.matrix.tsv .
    fi
    """
}


process EXTRACT_FASTA {
    publishDir "${params.outpath}/parsnp", mode: "copy"

    container "snads/xmfa-to-fasta@sha256:ef61b6d2c1a3ac675ecd102d64488f65f715745788deaf9fae2d7ab69c71c277"

    input:
        path(parsnp_xmfa)

    output:
        path("parsnp.fasta"), emit: parsnp_fasta

    script:
    """
    convert_xmfa_to_fasta.py --xmfa ${parsnp_xmfa} > parsnp.fasta
    sed -i.bak 's/.ref//g' parsnp.fasta
    rm parsnp.fasta.bak
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """

    stub:
    """
    if [[ ! -f parsnp.fasta ]]; then
        cp ${params.examplepath}/parsnp/parsnp.fasta .
    fi
    """
}


process INFER_RECOMBINATION_GUBBINS {
    publishDir "${params.outpath}/gubbins", mode: "copy"

    container = "snads/gubbins@sha256:391a980312096f96d976f4be668d4dea7dda13115db004a50e49762accc0ec62"

    input:
        path(extracted_fasta)
        path(parsnp_tree)

    output:
        path("gubbins.recombination_predictions.gff"), emit: recombination_positions
        path("gubbins.node_labelled.final_tree.tre"), emit: node_labelled_tree

    script:
    """
    # Log version info
    echo "INFER_RECOMBINATION_GUBBINS:" >> ${params.logpath}/versions.txt
    run_gubbins.py --version | sed 's/^/    /' >> ${params.logpath}/versions.txt

    run_gubbins.py --starting-tree ${parsnp_tree} --prefix gubbins ${extracted_fasta}
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """

    stub:
    """
    if [[ ! -f gubbins.recombination_predictions.gff || ! -f gubbins.node_labelled.final_tree.tre ]]; then
        cp ${params.examplepath}/gubbins/gubbins.recombination_predictions.gff .
        cp ${params.examplepath}/gubbins/gubbins.node_labelled.final_tree.tre .
    fi
    """
}


process INFER_RECOMBINATION_CFML {
    publishDir "${params.outpath}/clonalframeml", mode: "copy"

    container = "snads/clonalframeml@sha256:bc00db247195fdc6151793712a74cc9b272dc2c9f153bb0037415e387f15351e"

    input:
        path(extracted_fasta)
        path(parsnp_tree)

    output:
        path("clonalframeml.importation_status.txt"), emit: recombination_positions
        path("clonalframeml.labelled_tree.newick"), emit: node_labelled_tree

    script:
    """
    # Log version info
    echo "INFER_RECOMBINATION_CFML:" >> ${params.logpath}/versions.txt
    ClonalFrameML -version | sed 's/^/    /' >> ${params.logpath}/versions.txt

    # ClonalFrameML needs tree labels to not be surrounded by single quotes
    sed -i.bak "s/'//g" ${parsnp_tree}
    rm ${parsnp_tree}.bak
    ClonalFrameML ${parsnp_tree} ${extracted_fasta} clonalframeml
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """

    stub:
    """
    if [[ ! -f clonalframeml.importation_status.txt || ! -f clonalframeml.labelled_tree.newick ]]; then
        cp ${params.examplepath}/clonalframeml/clonalframeml.importation_status.txt .
        cp ${params.examplepath}/clonalframeml/clonalframeml.labelled_tree.newick .
    fi
    """
}


process MASK_RECOMBINATION {
    publishDir "${params.outpath}/${recombination_method}", mode: "copy"

    container = "snads/mask-recombination@sha256:0df4f5e26b2beeb2a180c2e4d75148cde55d4cc62585b5053d6606c6623d33e4"

    input:
        path(extracted_fasta)
        path(node_labelled_tree)
        path(recombination_positions)
        val recombination_method

    output:
        path("${recombination_method}_masked_recombination.fasta"), emit: masked_fasta

    script:
    """
    mask_recombination.py \
        --alignment ${extracted_fasta} \
        --format ${recombination_method} \
        --rec_positions ${recombination_positions} \
        --tree ${node_labelled_tree} \
        > ${recombination_method}_masked_recombination.fasta
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """

    stub:
    """
    if [[ ! -f ${recombination_method}_masked_recombination.fasta ]]; then
        cp ${params.examplepath}/${recombination_method}/${recombination_method}_masked_recombination.fasta .
    fi
    """
}


process REINFER_TREE {
    publishDir "${params.outpath}/${recombination_method}", mode: "copy"

    container "snads/parsnp@sha256:0dc552de1cf91fb30aa25eb09b4d9eef841abae989760c937b0006dacd165377"

    input:
        path(masked_fasta)
        val recombination_method

    output:
        path("${recombination_method}_masked_recombination.tree"), emit: reinferred_tree

    script:
    """
    # Log version info
    echo "REINFER_TREE:" >> ${params.logpath}/versions.txt
    parsnp --version | sed 's/^/    /' >> ${params.logpath}/versions.txt
    fasttree -expert &> tmp.txt; head -1 tmp.txt | sed 's/^/    /' >> ${params.logpath}/versions.txt; rm tmp.txt
    raxmlHPC-PTHREADS -v | sed 's/^/    /' >> ${params.logpath}/versions.txt; rm tmp.txt

    if [ "${params.reinferTreeProg}" = "fasttree" ]; then
        fasttree -nt ${masked_fasta} > ${recombination_method}_masked_recombination.tree
    elif [ "${params.reinferTreeProg}" = "raxml" ]; then
        raxmlHPC-PTHREADS -s ${masked_fasta} -m GTRGAMMA -w ${workDir} -n ${recombination_method}_masked_recombination -p 5280
        mv ${workDir}/RAxML_bestTree.${recombination_method}_masked_recombination ${recombination_method}_masked_recombination.tree
    fi
    if [ ! -s ${recombination_method}_masked_recombination.tree ]; then
        echo "Empty tree file. Did tree building program run out of RAM?" >> ${params.logpath}/stderr.nextflow.txt
        exit 1
    fi

    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """

    stub:
    """
    if [[ ! -f ${recombination_method}_masked_recombination.tree ]]; then
        cp ${params.examplepath}/${recombination_method}/${recombination_method}_masked_recombination.tree .
    fi
    """
}