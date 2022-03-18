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
    The minimal command for running the pipeline is:
    nextflow run main.nf
    A more typical command for running the pipeline is:
    nextflow run -profile singularity main.nf --inpath INPUT_DIR --outpath OUTPATH_DIR

    Input/output options:
      --inpath             Path to input data directory containing FastA assemblies. Recognized extensions are:  fa, fasta, fas, fna, fsa, fa.gz, fasta.gz, fas.gz, fna.gz, fsa.gz.
      --outpath            The output directory where the results will be saved.
    Analysis options:
      --curated-input      Whether or not input is a curated genome directory. If true, will assemue all genomes are similar enough to return sensible results. Options are: true (default), false.
      --recombination      Use a program to classify SNPs as due to recombination. Options are: gubbins, cfml, both.
      --reinfer-tree-prog  Program used to re-infer tree without SNPs classified as due to recombination. Options are: fasttree (default), raxml.
      --max-partition-size Max partition size (in bases, limits ParSNP memory usage). Note: results can change slightly depending on this value. Default is: 15000000.
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

// Handle command line input
if (!(params.recombination in ["cfml", "gubbins", "both", false])){
    System.err.println "\nERROR: --recombination was set as: ${params.recombination}"
    System.err.println "\nERROR: --recombination must be: cfml|gubbins|both"
    exit 1
}

if (!(params.curatedInput in ["true", "false", true, false])){
    System.err.println "\nERROR: --curated-input was set as: ${params.curatedInput}"
    System.err.println "\nERROR: --curated-input must be: true|false"
    exit 1
}

if (!(params.reinferTreeProg in ["fasttree", "raxml", false])){
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
    recombination:      ${params.recombination}
    curated-input:      ${params.curatedInput}
    max-partition-size: ${params.maxPartitionSize}
    reinfer-tree:       ${params.recombination ? params.reinferTreeProg : '-'}
    refpath:            ${params.refpath}
    =====================================
    """
    .stripIndent()

/*
==============================================================================
                 Import local custom modules and subworkflows                 
==============================================================================
*/
include {
    FIND_INFILES;
    INFILE_HANDLING;
    RUN_PARSNP;
    EXTRACT_SNPS;
    PAIRWISE_DISTANCES;
    DISTANCE_MATRIX;
    CLEANUP_FILES;
    EXTRACT_FASTA;
    INFER_RECOMBINATION_GUBBINS;
    INFER_RECOMBINATION_CFML;
    MASK_RECOMBINATION;
    REINFER_TREE
} from "./modules/assembly-snps.nf"

include { RUN_GUBBINS } from './subworkflows/run_gubbins.nf'
include { RUN_CFML } from './subworkflows/run_cfml.nf'

/*
==============================================================================
                   Import nf-core modules and subworkflows                    
==============================================================================
*/
            //
            // none
            //

/*
==============================================================================
                            Run the main workflow                             
==============================================================================
*/

workflow {

    // Generate a core-genome alignment, call SNPs, and build a tree
    inpath = Channel.fromPath(params.inpath, checkIfExists: true) // a filepath
    refpath = Channel.from(params.refpath) // a string (filepath or 'largest')
    outpath = Channel.fromPath(params.outpath)

    FIND_INFILES(
        inpath
    )

    INFILE_HANDLING(
        FIND_INFILES.out.find_infiles_success,
        inpath,
        refpath
    )

    RUN_PARSNP(
        INFILE_HANDLING.out.refpath,
        INFILE_HANDLING.out.tmppath
    )

    EXTRACT_SNPS(
        RUN_PARSNP.out.parsnp_ggr
    )

    PAIRWISE_DISTANCES(
        EXTRACT_SNPS.out.snps_file
    )

    DISTANCE_MATRIX(
        PAIRWISE_DISTANCES.out.snp_distances
    )

    CLEANUP_FILES(
        PAIRWISE_DISTANCES.out.snp_distances,
        outpath
    )

    // Optionally infer SNPs due to recombination, mask them, and build a new tree
    if (params.recombination != false) {
        EXTRACT_FASTA(
            RUN_PARSNP.out.parsnp_xmfa
        )
    }

    if (params.recombination == 'gubbins') {
        RUN_GUBBINS(
            EXTRACT_FASTA.out.parsnp_fasta,
            RUN_PARSNP.out.parsnp_tree
        )
    } else if (params.recombination == 'cfml') {
        RUN_CFML(
            EXTRACT_FASTA.out.parsnp_fasta,
            RUN_PARSNP.out.parsnp_tree
        )
    } else if (params.recombination == 'both') {
        RUN_GUBBINS(
            EXTRACT_FASTA.out.parsnp_fasta,
            RUN_PARSNP.out.parsnp_tree
        )
        RUN_CFML(
            EXTRACT_FASTA.out.parsnp_fasta,
            RUN_PARSNP.out.parsnp_tree
        )
    }
}

/*
==============================================================================
                        Completion e-mail and summary                         
==============================================================================
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
