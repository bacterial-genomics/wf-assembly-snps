#!/usr/bin/env nextflow


/*
==============================================================================
                              wf-assembly-snps                              
==============================================================================
usage: nextflow run ./wf-assembly-snps/main.nf [-help]
----------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    =========================================
     wf-assembly-snps v${version}
    =========================================

    Usage:
    The minimal command for running the workflow is:
    nextflow run main.nf
    To run the workflow on a small set of test data:
    nextflow run main.nf -profile test,<docker|singularity> --outpath results
    A typical command for running the pipeline is:
    nextflow run -profile <docker|singularity> main.nf --inpath <input directory> --outpath <directory for results>

    Input/output options:
      --inpath             Path to input data directory containing FastA assemblies. Recognized extensions are:  fa, fasta, fas, fna, fsa, fa.gz, fasta.gz, fas.gz, fna.gz, fsa.gz.
      --outpath            The output directory where the results will be saved.
    Analysis options:
      --curated_input      Whether or not input is a curated genome directory. If true, will assemue all genomes are similar enough to return sensible results. Options are: true (default), false.
      --recombination      Use a program to classify SNPs as due to recombination. Options are: gubbins, cfml, both.
      --tree_method        Program used to infer trees (in ParSNP and optionally again after masking positions due to recombination). Options are: fasttree (default), raxml.
      --max_partition_size Max partition size (in bases, limits ParSNP memory usage). Note: results can change slightly depending on this value. Default is: 15000000.
      --bigdata            Whether or not to use more compute resources. Options are true, false (default).
      --max_memory         Specify memory limit on your machine/infrastructure, e.g. '128.GB'. Useful to ensure workflow doesn't request too many resources.
      --max_time           Specify time limit for each process, e.g. '240.h'. Useful to ensure workflow doesn't request too many resources.
      --max_cpus           Specify CPU limit on your machine/infrastructure, e.g. 16. Useful to ensure workflow doesn't request too many resources.
      --min_ggr_size       Minimum filesize to verify that ParSNP produced a .ggr file. Default is: '100k'.
      --min_xmfa_size      Minimum filesize to verify that ParSNP produced a .xmfa file. Default is: '100k'.
    Profile options:
      -profile singularity Use Singularity images to run the workflow. Will pull and convert Docker images from Dockerhub if not locally available.
      -profile docker      Use Docker images to run the workflow. Will pull images from Dockerhub if not locally available.
      -profile conda       TODO: this is not implemented yet.
    Other options:
      -resume              Re-start a workflow using cached results. May not behave as expected with containerization profiles docker or singularity.
      -stub                Use example output files for any process with an uncommented stub block. For debugging/testing purposes.
      -name                Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
    """.stripIndent()
}

version = "1.0.0"
nextflow.enable.dsl=2

if (params.help) {
    helpMessage()
    exit 0
}

if (params.version){
    println "VERSION: $version"
    exit 0
}

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

if (!(params.recombination in ["cfml", "gubbins", "both", false])){
    System.err.println "\nERROR: --recombination was set as: ${params.recombination}"
    System.err.println "\nERROR: --recombination must be: cfml|gubbins|both"
    exit 1
}

if (!(params.curated_input in ["true", "false", true, false])){
    System.err.println "\nERROR: --curated-input was set as: ${params.curatedInput}"
    System.err.println "\nERROR: --curated-input must be: true|false"
    exit 1
}

if (!(params.tree_method in ["fasttree", "raxml", false])){
    System.err.println "\nERROR: --reinfer-tree-prog was set as: ${params.reinferTreeProg}"
    System.err.println "\nERROR: --reinfer-tree-prog must be: fasttree|raxml"
    exit 1
}

File inpathFileObj = new File(params.inpath)
if (!inpathFileObj.exists()){
    System.err.println "ERROR: $params.inpath doesn't exist"
    exit 1
}

File outpathFileObj = new File(params.outpath)
if (outpathFileObj.exists()){
    // Per the config file, outpath stores log & trace files so it is created before this point
    // Check that outpath only contains a trace file created this hour
    dayAndHour = new java.util.Date().format('yyyy-MM-dd HH')
    outFiles = outpathFileObj.list()
    if (!(outFiles[0] ==~ /trace.($dayAndHour):\d\d:\d\d.txt/ && outFiles.size() == 1)) {
        // If it contains an older trace file or other files, warn the user
        System.out.println "WARNING: $params.outpath already exists. Output files will be overwritten."
    }
} else {
    outpathFileObj.mkdirs()
}

File logpathFileObj = new File(params.logpath)
if (logpathFileObj.exists()){
    System.out.println "WARNING: $params.logpath already exists. Log files will be overwritten."
} else {
    logpathFileObj.mkdirs()
}

// Print parameters used
log.info """
    =====================================
    wf-assembly-snps $version
    =====================================
    inpath:             ${params.inpath}
    outpath:            ${params.outpath}
    logpath:            ${params.logpath}
    workDir:            ${workflow.workDir}
    recombination:      ${params.recombination}
    curated_input:      ${params.curated_input}
    max_partition_size: ${params.max_partition_size}
    tree_method:        ${params.tree_method}
    refpath:            ${params.refpath}
    =====================================
    """
    .stripIndent()

/*
========================================================================================
                 Import local custom modules and subworkflows                 
========================================================================================
*/
include { INFILE_HANDLING } from "./modules/local/infile_handling"
include { PARSNP } from "./modules/local/parsnp"
include { EXTRACT_SNPS } from "./modules/local/extract_snps"
include { PAIRWISE_DISTANCES } from "./modules/local/pairwise_distances"
include { DISTANCE_MATRIX } from "./modules/local/distance_matrix"
include { OUTFILE_CLEANUP } from "./modules/local/outfile_cleanup"
include { EXTRACT_FASTA } from "./modules/local/extract_fasta"

include { RECOMBINATION_GUBBINS } from './subworkflows/recombination_gubbins'
include { RECOMBINATION_CFML } from './subworkflows/recombination_cfml'

/*
========================================================================================
                   Import nf-core modules and subworkflows                    
========================================================================================
*/

// None

/*
========================================================================================
                            Run the main workflow                             
========================================================================================
*/

workflow {

    // SETUP: Define input & output channels
    ch_inpath = Channel.fromPath(params.inpath, checkIfExists: true)
    ch_outpath = Channel.fromPath(params.outpath)
    if (params.refpath) {
        ch_refpath = Channel.fromPath(params.refpath, checkIfExists: true)
    } else {
        ch_refpath = Channel.fromPath('largest')  // dummy filepath hack (https://github.com/nextflow-io/nextflow/issues/1532)
    }
    ch_versions = Channel.empty()

    // PROCESS: Read in input directory path, validate and stage input files
    INFILE_HANDLING(
        ch_inpath,
        ch_refpath
    )
    ch_versions = ch_versions.mix(INFILE_HANDLING.out.versions)

    // PROCESS: Run ParSNP to generate core genome alignment and phylogeny
    PARSNP(
        INFILE_HANDLING.out.reference,
        INFILE_HANDLING.out.tmpdir
    )
    ch_versions = ch_versions.mix(PARSNP.out.versions)

    // PROCESS: Get SNP alignment from ParSNP .ggr file
    EXTRACT_SNPS(
        PARSNP.out.ggr
    )
    ch_versions = ch_versions.mix(EXTRACT_SNPS.out.versions)

    // PROCESS: Calculate pairwise genome distances
    PAIRWISE_DISTANCES(
        EXTRACT_SNPS.out.snps_file
    )
    ch_versions = ch_versions.mix(PAIRWISE_DISTANCES.out.versions)

    // PROCESS: Reformat pairwise genome distances into matrix
    DISTANCE_MATRIX(
        PAIRWISE_DISTANCES.out.snp_distances
    )
    ch_versions = ch_versions.mix(DISTANCE_MATRIX.out.versions)

    // PROCESS: Gzip compress SNP alignment
    OUTFILE_CLEANUP(
        PAIRWISE_DISTANCES.out.snp_distances,  // process waits until this file is produced
        ch_outpath
    )
    ch_versions = ch_versions.mix(OUTFILE_CLEANUP.out.versions)

    // PROCESS: Get core genome alignment from ParSNP .xmfa file
    if (params.recombination != false) {
        EXTRACT_FASTA(
            PARSNP.out.xmfa
        )
        ch_versions = ch_versions.mix(EXTRACT_FASTA.out.versions)
    }

    // WORKFFLOW: Infer SNPs due to recombination, mask them, re-infer phylogeny
    if (params.recombination == 'gubbins') {
        RECOMBINATION_GUBBINS(
            EXTRACT_FASTA.out.parsnp_fasta,
            PARSNP.out.tree
        )
    } else if (params.recombination == 'cfml') {
        RECOMBINATION_CFML(
            EXTRACT_FASTA.out.parsnp_fasta,
            PARSNP.out.tree
        )
        ch_versions = ch_versions.mix(RECOMBINATION_CFML.out.cfml_versions)
        ch_versions = ch_versions.mix(RECOMBINATION_CFML.out.mask_recombination_versions)
        ch_versions = ch_versions.mix(RECOMBINATION_CFML.out.infer_tree_versions)
    } else if (params.recombination == 'both') {
        RECOMBINATION_GUBBINS(
            EXTRACT_FASTA.out.parsnp_fasta,
            PARSNP.out.tree
        )
        RECOMBINATION_CFML(
            EXTRACT_FASTA.out.parsnp_fasta,
            PARSNP.out.tree
        )
    }

    // PATTERN: Collate method version information
    ch_versions.collectFile(name: 'software_versions.yml', storeDir: params.outpath)

}

/*
========================================================================================
                        Completion e-mail and summary                         
========================================================================================
*/

workflow.onComplete {
    log.info """
                |=====================================
                |Pipeline Execution Summary
                |=====================================
                |Workflow Version : ${version}
                |Nextflow Version : ${nextflow.version}
                |Command Line     : ${workflow.commandLine}
                |Resumed          : ${workflow.resume}
                |Completed At     : ${workflow.complete}
                |Duration         : ${workflow.duration}
                |Success          : ${workflow.success}
                |Exit Code        : ${workflow.exitStatus}
                |Launch Dir       : ${workflow.launchDir}
                |=====================================
             """.stripMargin()
}

workflow.onError {
    def err_msg = """
                     |=====================================
                     |Error summary
                     |=====================================
                     |Completed at : ${workflow.complete}
                     |exit status  : ${workflow.exitStatus}
                     |workDir      : ${workflow.workDir}
                     |Error Report :
                     |${workflow.errorReport ?: '-'}
                     |=====================================
                  """.stripMargin()
    log.info err_msg

/*    sendMail(
        to: "${USER}@cdc.gov",
        subject: 'workflow error',
        body: err_msg,
        charset: UTF-8
        // attach: '/path/stderr.log.txt'
    )
*/
}
