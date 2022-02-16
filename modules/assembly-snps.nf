nextflow.enable.dsl = 2


process FIND_INFILES {

    input:
        path(input_path)

    output:
        path("find_infiles.success.txt"), emit: found_infiles

    script:
    """
    find_infiles.sh ${input_path}
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """

//     stub:
//     """
//     touch find_infiles.success.txt
//     """

}


process INFILE_HANDLING {
    
    input:
        path(found_infiles)
        path(input_dir_path)
        path(output_dir_path)
        path(reference_file_path)

    output:
        path("${output_dir_path}/.ref/*"), emit: ref_path
        path("${output_dir_path}/.tmp"), emit: tmp_path
        path("${output_dir_path}/.tmp/*")
        path("infile_handling.success.txt"), emit: handled_infiles

    script:
    """
    infile_handling.sh ${input_dir_path} ${output_dir_path} ${reference_file_path}
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """
    
    stub:
    """
    if [[ ! -d  ${output_dir_path}/.tmp || ${output_dir_path}/.ref ]]; then
        cp -r ${params.examplepath}/.tmp ${output_dir_path}
        cp -r ${params.examplepath}/.ref ${output_dir_path}
    fi
    touch infile_handling.success.txt
    """

}


process RUN_PARSNP {

    cpus 2
    params.enable_conda_yml ? "$baseDir/conda/linux/parsnp.yml" : null
    //conda 'bioconda::parsnp=1.1.3'
    container "snads/parsnp:1.5.6"
    // @sha256:f43ffe7ed111c9721891950d25160d94cec3d8dbdaf4f3afc16ce91d184ed34f
    // process.container = "dockerhub_user/image_name:image_tag"
    //container "hub.docker.com/layers/snads/parsnp/1.5.6"

    input:
        path(handled_infiles)
        path(ref_path)
        path(tmp_path)
        path(output_dir_path)

    output:
        path("${output_dir_path}/parsnp.ggr")
        path("${output_dir_path}/parsnp.xmfa"), emit: parsnp_xmfa
        path("${output_dir_path}/parsnp.tree"), emit: parsnp_tree
        path("run_parsnp.success.txt"), emit: ran_parsnp

    script:
    """
    run_parsnp.sh ${tmp_path} ${output_dir_path} ${ref_path} ${task.cpus}
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """

    stub:
    """
    if [[ ! -d  ${output_dir_path}/parsnp.ggr || ${output_dir_path}/parsnp.tree || ${output_dir_path}/parsnp.xmfa ]]; then
        cp ${params.examplepath}/parsnp.ggr  ${output_dir_path}
        cp ${params.examplepath}/parsnp.xmfa  ${output_dir_path}
        cp ${params.examplepath}/parsnp.tree ${output_dir_path}
    fi
    touch run_parsnp.success.txt
    """


}

process EXTRACT_SNPS {


    params.enable_conda_yml ? "$baseDir/conda/linux/harvesttools.yml" : null
    // conda 'bioconda::harvesttools=1.2'
    // container = "$baseDir/assets/parsnp_1.5.6.sif"  // TODO: replace placeholder with option to run with singularity
    container "snads/parsnp:1.5.6"

    input:
        path(ran_parsnp)
        path(output_dir_path)

    output:
        path("${output_dir_path}/SNPs.fa")
        path("extract_snps.success.txt"), emit: extracted_snps

    script:
    """
    extract_snps.sh ${output_dir_path}
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """

    stub:
    """
    if [[ ! -f ${output_dir_path}/SNPs.fa ]]; then
        cp ${params.examplepath}/SNPs.fa  ${output_dir_path}
    fi
    touch extract_snps.success.txt
    """

}

process PAIRWISE_DISTANCES {

    params.enable_conda_yml ? "$baseDir/conda/linux/NEEDS-NEWFILE.yml" : null
    // conda 'bioconda::FIXME'
    // container = "$baseDir/assets/parsnp_1.5.6.sif"  // TODO: replace placeholder with option to run with singularity
    container "snads/hamming-dist:1.0"
    cpus 2

    input:
        path(extracted_snps)
        path(output_dir_path)

    output:
        path("${output_dir_path}/SNP-distances.pairs.tsv")
        path("pairwise_distances.success.txt"), emit: calculated_snp_distances

    script:
    """
    pairwise_distances.sh ${output_dir_path} ${task.cpus}
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """

    stub:
    """
    if [[ ! -f ${output_dir_path}/SNP-distances.pairs.tsv ]]; then
        cp ${params.examplepath}/SNP-distances.pairs.tsv  ${output_dir_path}
    fi
    touch pairwise_distances.success.txt
    """

}


process DISTANCE_MATRIX {

    params.enable_conda_yml ? "$baseDir/conda/linux/python3.yml" : null
    // conda 'conda-forge::python=3.10.1'
    // container = "$baseDir/assets/python_3.sif"  // TODO: replace placeholder with option to run with singularity
    container "snads/hamming-dist:1.0"

    input:
        path(calculated_snp_distances)
        path(output_dir_path)

    output:
        path("${output_dir_path}/SNP-distances.matrix.tsv")
        path("distance_matrix.success.txt"), emit: made_snp_matrix

    script:
    """
    distance_matrix.sh ${output_dir_path}
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """

    stub:
    """
    if [[ ! -f ${output_dir_path}/SNP-distances.matrix.tsv ]]; then
        cp ${params.examplepath}/SNP-distances.matrix.tsv  ${output_dir_path}
    fi
    touch distance_matrix.success.txt
    """
}

process EXTRACT_FASTA {
    container "snads/xmfa-to-fasta:2.0"

    input:
        path(parsnp_xmfa)
        path(output_dir_path)

    output:
        path("${output_dir_path}/parsnp.fasta"), emit: parsnp_fasta
        path("extract_fasta.success.txt"), emit: extracted_fasta

    script:
    """
    if convert_xmfa_to_fasta.py --xmfa ${parsnp_xmfa} > ${output_dir_path}/parsnp.fasta; then
        sed -i.bak 's/.ref//g' ${output_dir_path}/parsnp.fasta
        rm ${output_dir_path}/parsnp.fasta.bak
        touch "extract_fasta.success.txt"
    fi
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """

    stub:
    """
    if [[ ! -f ${output_dir_path}/parsnp.fasta ]]; then
        cp ${params.examplepath}/parsnp.fasta  ${output_dir_path}
    fi
    touch extract_fasta.success.txt
    """
}

process INFER_RECOMBINATION_GUBBINS {
    container = "snads/gubbins:3.1.4"

    input:
        path(extracted_fasta)
        path(parsnp_tree)
        path(output_dir_path)

    output:
        path("${output_dir_path}/parsnp.recombination_predictions.gff"), emit: recombination_positions
        path("infer_recombination.success.txt"), emit: inferred_recombination

    script:
    """
    if run_gubbins.py --starting-tree ${parsnp_tree} --prefix ${output_dir_path}/parsnp ${extracted_fasta}; then
        touch infer_recombination.success.txt
    fi
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """

    stub:
    """
    if [[ ! -f ${output_dir_path}/parsnp.recombination_predictions.gff ]]; then
        cp ${params.examplepath}/parsnp.recombination_predictions.gff  ${output_dir_path}
    fi
    touch infer_recombination.success.txt
    """
}

process INFER_RECOMBINATION_CFML {
    container = "snads/clonalframeml:1.12"

    input:
        path(extracted_fasta)
        path(parsnp_tree)
        path(output_dir_path)

    output:
        path("${output_dir_path}/parsnp.importation_status.txt"), emit: recombination_positions
        path("infer_recombination.success.txt"), emit: inferred_recombination

    script:
    """
    # ClonalFrameML needs tree labels to not be surrounded by single quotes
    sed -i.bak "s/'//g" ${parsnp_tree}
    rm ${parsnp_tree}.bak
    if ClonalFrameML ${parsnp_tree} ${extracted_fasta} ${output_dir_path}/parsnp; then
        touch infer_recombination.success.txt
    fi
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """

    stub:
    """
    if [[ ! -f ${output_dir_path}/parsnp.importation_status.txt ]]; then
        cp ${params.examplepath}/parsnp.importation_status.txt  ${output_dir_path}
    fi
    touch infer_recombination.success.txt
    """
}

process MASK_RECOMBINATION {
    container = "snads/mask-recombination:1.0"

    input:
        path(extracted_fasta)
        path(parsnp_tree)
        path(recombination_positions)
        val format
        path(output_dir_path)

    output:
        path("${output_dir_path}/parsnp_${format}.fasta"), emit: masked_fasta
        path("mask_recombination_${format}.success.txt"), emit: masked_recombination

    script:
    """
    if mask_recombination.py --alignment ${extracted_fasta} --format ${format} --rec_positions ${recombination_positions} --tree ${parsnp_tree} > ${output_dir_path}/parsnp_${format}.fasta; then
        touch mask_recombination_${format}.success.txt
    fi
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """

    stub:
    """
    if [[ ! -f ${output_dir_path}/parsnp_${format}.fasta ]]; then
        cp ${params.examplepath}/parsnp_${format}.fasta  ${output_dir_path}
    fi
    touch mask_recombination_${format}.success.txt
    """
}

process REINFER_TREE {
    container "snads/parsnp:1.5.6"

    input:
        path(masked_fasta)
        val recombination_method
        path(output_dir_path)

    output:
        path("${output_dir_path}/parsnp_${recombination_method}.tree")
        path("reinfer_tree_${recombination_method}.success.txt"), emit: reinferred_tree

    script:
    """
    # With RAxML
    # raxmlHPC-PTHREADS -s ${masked_fasta} -m GTRGAMMA -w ${launchDir}/${output_dir_path} -n parsnp_${recombination_method} -p 5280
    # mv ${launchDir}/${output_dir_path}/RAxML_result.parsnp_${recombination_method} ${output_dir_path}/parsnp_${recombination_method}.tree
    # With FastTree
    fasttree -nt ${masked_fasta} > ${output_dir_path}/parsnp_${recombination_method}.tree
    if [ -f "${output_dir_path}/parsnp_${recombination_method}.tree" ] && \
        [ -s "${output_dir_path}/parsnp_${recombination_method}.tree" ]; then
        touch reinfer_tree_${recombination_method}.success.txt
    fi
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """

//     stub:
//     """
//     if [[ ! -f ${output_dir_path}/parsnp_masked_recombination.tree ]]; then
//         cp ${params.examplepath}/parsnp_masked_recombination.tree  ${output_dir_path}
//     fi
//     touch reinfer_tree_${recombination_method}.success.txt
//     """
}