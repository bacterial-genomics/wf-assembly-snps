#!/usr/bin/env nextflow


/*
==============================================================================
                              wf-assembly-snps                              
==============================================================================
usage: nextflow run ./wf-assembly-snps/main.nf [-help]
----------------------------------------------------------------------------
*/

version = "1.0.0"
nextflow.enable.dsl=2

// Default parameters (can be overwritten by command line inputs)
params.inpath = new File("${launchDir}").getCanonicalPath()
params.outpath = new File("${launchDir}").getCanonicalPath()
params.logpath = new File("${params.outpath}/.log").getCanonicalPath()
params.refpath = 'largest'
params.examplepath = new File("${launchDir}/example_output").getCanonicalPath() // TODO: remove this when development is finished
params.help = false
params.version = false
params.recombination = false
params.enable_conda_yml = false

// Handle command line input
if (!(params.recombination in ["cfml", "gubbins", "both", false])){
    System.err.println "\nERROR: --recombination was set as: ${params.recombination}"
    System.err.println "\nERROR: --recombination must be: cfml|gubbins|both"
    exit 1
}

if (params.help){
    println "USAGE: put usage info here"
    exit 0
}

if (params.version){
    println "VERSION: $version"
    exit 0
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
    inpath:         ${params.inpath}
    outpath:        ${params.outpath}
    logpath:        ${params.logpath}
    recombination:  ${params.recombination}
    refpath:        ${params.refpath}
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
