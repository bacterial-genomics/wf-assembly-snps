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

    stub:
    """
    touch find_infiles.success.txt
    """

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
    touch infile_handling.success.txt
    """

}


process RUN_PARSNP {

    cpus 2
    conda (params.enable_conda_yml ? "$baseDir/conda/linux/parsnp.yml" : 'bioconda::parsnp=1.5.6')
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
        path("${output_dir_path}/parsnp.tree")
        path("run_parsnp.success.txt"), emit: ran_parsnp

    script:
    """
    run_parsnp.sh ${tmp_path} ${output_dir_path} ${ref_path} ${task.cpus}
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """

    stub:
    """
    touch run_parsnp.success.txt
    """


}

process EXTRACT_SNPS {


    conda (params.enable_conda_yml ? "$baseDir/conda/linux/harvesttools.yml" : 'bioconda::harvesttools=1.2')
    // conda (params.enable_conda_yml ? 'bioconda::harvesttools=1.2' : null)
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
    touch extract_snps.success.txt
    """

}

process PAIRWISE_DISTANCES {

    conda (params.enable_conda_yml ? "$baseDir/conda/linux/bioperl.yml" : 'bioconda::perl-bioperl-core=1.007002')
    // conda (params.enable_conda_yml ? 'bioconda::perl-bioperl-core=1.007002' : null)
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
    touch pairwise_distances.success.txt
    """

}


process DISTANCE_MATRIX {

    conda (params.enable_conda_yml ? "$baseDir/conda/linux/python3.yml" : 'conda-forge::python=3.10.1')
    // conda (params.enable_conda_yml ? 'conda-forge::python=3.10.1' : null)
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
    touch distance_matrix.success.txt
    """
}


process RECOMBINATION {

    cpus 2
    container = "$baseDir/assets/filename.sif"

    input:
        path(extracted_snps)
        path(output_dir_path)

    output:
        path("${output_dir_path}/mutational-only.recombination-free.tre")
        path("recombination.success.txt"), emit: inferred_recombination

    script:
    """
    recombination.sh ${output_dir_path} ${task.cpus}
    cat .command.out >> ${params.logpath}/stdout.nextflow.txt
    cat .command.err >> ${params.logpath}/stderr.nextflow.txt
    """

    stub:
    """
    touch recombination.success.txt
    """
}
