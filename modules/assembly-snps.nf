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
        path(outpath)
        val refpath

    output:
        path("${outpath}/.ref/*"), emit: refpath
        path("${outpath}/.tmp"), emit: tmppath
        path("${outpath}/.tmp/*")
        path("infile_handling.success.txt"), emit: infile_handling_success

    script:
    """
    infile_handling.sh ${inpath} ${outpath} ${refpath}
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """
    
    stub:
    """
    if [[ ! -d  ${outpath}/.tmp || ${outpath}/.ref ]]; then
        cp -r ${params.examplepath}/.tmp ${outpath}
        cp -r ${params.examplepath}/.ref ${outpath}
    fi
    touch infile_handling.success.txt
    """

}


process RUN_PARSNP {

    params.enable_conda_yml ? "$baseDir/conda/linux/parsnp.yml" : null
    //conda 'bioconda::parsnp=1.1.3'
    container "snads/parsnp:1.5.6"

    input:
        path(infile_handling_success)
        path(refpath)
        path(tmppath)
        path(outpath)

    output:
        path("${outpath}/parsnp.ggr")
        path("${outpath}/parsnp.xmfa"), emit: parsnp_xmfa
        path("${outpath}/parsnp.tree"), emit: parsnp_tree
        path("run_parsnp.success.txt"), emit: run_parsnp_success

    script:
    """
    run_parsnp.sh ${tmppath} ${outpath} ${refpath} ${task.cpus}
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """

    stub:
    """
    if [[ ! -d  ${outpath}/parsnp.ggr || ${outpath}/parsnp.tree || ${outpath}/parsnp.xmfa ]]; then
        cp ${params.examplepath}/parsnp.ggr  ${outpath}
        cp ${params.examplepath}/parsnp.xmfa  ${outpath}
        cp ${params.examplepath}/parsnp.tree ${outpath}
    fi
    touch run_parsnp.success.txt
    """

}


process EXTRACT_SNPS {

    params.enable_conda_yml ? "$baseDir/conda/linux/harvesttools.yml" : null
    // conda 'bioconda::harvesttools=1.2'
    container "snads/parsnp:1.5.6"

    input:
        path(run_parsnp_success)
        path(outpath)

    output:
        path("${outpath}/SNPs.fa")
        path("extract_snps.success.txt"), emit: extract_snps_success

    script:
    """
    extract_snps.sh ${outpath}
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """

    stub:
    """
    if [[ ! -f ${outpath}/SNPs.fa ]]; then
        cp ${params.examplepath}/SNPs.fa  ${outpath}
    fi
    touch extract_snps.success.txt
    """

}


process PAIRWISE_DISTANCES {

    params.enable_conda_yml ? "$baseDir/conda/linux/NEEDS-NEWFILE.yml" : null
    // conda 'bioconda::FIXME'
    container "snads/hamming-dist:1.0"

    input:
        path(extract_snps_success)
        path(outpath)

    output:
        path("${outpath}/SNP-distances.pairs.tsv")
        path("pairwise_distances.success.txt"), emit: pairwise_distances_success

    script:
    """
    pairwise_distances.sh ${outpath} ${task.cpus}
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """

    stub:
    """
    if [[ ! -f ${outpath}/SNP-distances.pairs.tsv ]]; then
        cp ${params.examplepath}/SNP-distances.pairs.tsv  ${outpath}
    fi
    touch pairwise_distances.success.txt
    """

}


process DISTANCE_MATRIX {

    params.enable_conda_yml ? "$baseDir/conda/linux/python3.yml" : null
    // conda 'conda-forge::python=3.10.1'
    container "snads/hamming-dist:1.0"

    input:
        path(pairwise_distances_success)
        path(outpath)

    output:
        path("${outpath}/SNP-distances.matrix.tsv")
        path("distance_matrix.success.txt"), emit: distance_matrix_success

    script:
    """
    distance_matrix.sh ${outpath}
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """

    stub:
    """
    if [[ ! -f ${outpath}/SNP-distances.matrix.tsv ]]; then
        cp ${params.examplepath}/SNP-distances.matrix.tsv  ${outpath}
    fi
    touch distance_matrix.success.txt
    """
}


process EXTRACT_FASTA {

    container "snads/xmfa-to-fasta:2.0"

    input:
        path(run_parsnp_success)
        path(parsnp_xmfa)
        path(outpath)

    output:
        path("${outpath}/parsnp.fasta"), emit: parsnp_fasta
        path("extract_fasta.success.txt"), emit: extract_fasta_success

    script:
    """
    if convert_xmfa_to_fasta.py --xmfa ${parsnp_xmfa} > ${outpath}/parsnp.fasta; then
        sed -i.bak 's/.ref//g' ${outpath}/parsnp.fasta
        rm ${outpath}/parsnp.fasta.bak
        touch "extract_fasta.success.txt"
    fi
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """

    stub:
    """
    if [[ ! -f ${outpath}/parsnp.fasta ]]; then
        cp ${params.examplepath}/parsnp.fasta  ${outpath}
    fi
    touch extract_fasta.success.txt
    """
}


process INFER_RECOMBINATION_GUBBINS {

    container = "snads/gubbins:3.1.4"

    input:
        path(extracted_fasta)
        path(parsnp_tree)
        path(outpath)

    output:
        path("${outpath}/parsnp.recombination_predictions.gff"), emit: recombination_positions
        path("${outpath}/parsnp.node_labelled.final_tree.tre"), emit: node_labelled_tree
        path("infer_recombination.success.txt"), emit: infer_recombination_success

    script:
    """
    if run_gubbins.py --starting-tree ${parsnp_tree} --prefix ${outpath}/parsnp ${extracted_fasta}; then
        touch infer_recombination.success.txt
    fi
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """

    stub:
    """
    if [[ ! -f ${outpath}/parsnp.recombination_predictions.gff || ! -f ${outpath}/parsnp.node_labelled.final_tree.tre ]]; then
        cp ${params.examplepath}/parsnp.recombination_predictions.gff  ${outpath}
        cp ${params.examplepath}/parsnp.node_labelled.final_tree.tre  ${outpath}
    fi
    touch infer_recombination.success.txt
    """
}


process INFER_RECOMBINATION_CFML {

    container = "snads/clonalframeml:1.12"

    input:
        path(extracted_fasta)
        path(parsnp_tree)
        path(outpath)

    output:
        path("${outpath}/parsnp.importation_status.txt"), emit: recombination_positions
        path("${outpath}/parsnp.labelled_tree.newick"), emit: node_labelled_tree
        path("infer_recombination.success.txt"), emit: infer_recombination_success

    script:
    """
    # ClonalFrameML needs tree labels to not be surrounded by single quotes
    sed -i.bak "s/'//g" ${parsnp_tree}
    rm ${parsnp_tree}.bak
    if ClonalFrameML ${parsnp_tree} ${extracted_fasta} ${outpath}/parsnp; then
        touch infer_recombination.success.txt
    fi
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """

    stub:
    """
    if [[ ! -f ${outpath}/parsnp.importation_status.txt || ! -f ${outpath}/parsnp.labelled_tree.newick ]]; then
        cp ${params.examplepath}/parsnp.importation_status.txt  ${outpath}
        cp ${params.examplepath}/parsnp.labelled_tree.newick ${outpath}
    fi
    touch infer_recombination.success.txt
    """
}


process MASK_RECOMBINATION {

    container = "snads/mask-recombination:1.0"

    input:
        path(infer_recombination_success)
        path(extracted_fasta)
        path(node_labelled_tree)
        path(recombination_positions)
        val format
        path(outpath)

    output:
        path("${outpath}/parsnp_${format}.fasta"), emit: masked_fasta
        path("mask_recombination_${format}.success.txt"), emit: mask_recombination_success

    script:
    """
    if mask_recombination.py \
        --alignment ${extracted_fasta} \
        --format ${format} \
        --rec_positions ${recombination_positions} \
        --tree ${node_labelled_tree} \
        > ${outpath}/parsnp_${format}.fasta; then
        touch mask_recombination_${format}.success.txt
    fi
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """

    stub:
    """
    if [[ ! -f ${outpath}/parsnp_${format}.fasta ]]; then
        cp ${params.examplepath}/parsnp_${format}.fasta  ${outpath}
    fi
    touch mask_recombination_${format}.success.txt
    """
}


process REINFER_TREE {

    container "snads/parsnp:1.5.6"

    input:
        path(mask_recombination_success)
        path(masked_fasta)
        val recombination_method
        path(outpath)

    output:
        path("${outpath}/parsnp_${recombination_method}.tree")
        path("reinfer_tree_${recombination_method}.success.txt"), emit: reinfer_tree_success

    script:
    """
    # With RAxML
    raxmlHPC-PTHREADS -s ${masked_fasta} -m GTRGAMMA -w ${launchDir}/${outpath} -n parsnp_${recombination_method} -p 5280
    mv ${launchDir}/${outpath}/RAxML_result.parsnp_${recombination_method} ${outpath}/parsnp_${recombination_method}.tree
    # With FastTree
    # fasttree -nt ${masked_fasta} > ${outpath}/parsnp_${recombination_method}.tree
    if [ -f "${outpath}/parsnp_${recombination_method}.tree" ] && \
        [ -s "${outpath}/parsnp_${recombination_method}.tree" ]; then
        touch reinfer_tree_${recombination_method}.success.txt
    fi
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """

}