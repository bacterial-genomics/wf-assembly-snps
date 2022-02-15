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

params.help = false
if (params.help){
    println "USAGE: put usage info here"
    exit 0
}
params.version = false
if (params.version){
    println "VERSION: $version"
    exit 0
}

// Default parameters
params.inpath = new File("${launchDir}").getCanonicalPath()
params.outpath = new File("${launchDir}").getCanonicalPath()
params.logpath = new File("${params.outpath}/.log").getCanonicalPath()
params.refpath = new File("${launchDir}/INPUT_DIR/16-090.fna.gz").getCanonicalPath()
params.examplepath = new File("${launchDir}/example_output").getCanonicalPath()
//params.refpath = null
params.recombination = false

params.enable_conda_yml = false


// Checks on recombination parameter
if (!(params.recombination in ["cfml", "gubbins", "both", false])){
    System.err.println "\nERROR: --recombination was set as: ${params.recombination}"
    System.err.println "\nERROR: --recombination must be: cfml|gubbins|both"
    exit 1
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
         reference:      ${params.refpath}
         =====================================
         """
         .stripIndent()

// Path handling
File inpathFileObj = new File(params.inpath)
if (!inpathFileObj.exists()){
    System.err.println "ERROR: $params.inpath doesn't exist"
    exit 1
}
File outpathFileObj = new File(params.outpath)
if (outpathFileObj.exists()){
    System.out.println "WARNING: $params.outpath already exists. Output files will be overwritten."
} else {
    outpathFileObj.mkdirs()
}
File logpathFileObj = new File(params.logpath)
if (logpathFileObj.exists()){
    System.out.println "WARNING: $params.logpath already exists. Log files will be overwritten."
} else {
    logpathFileObj.mkdirs()
}


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
} from "./modules.nf"


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
    if (params.refpath) {
        ref_ch = Channel.fromPath(params.refpath, checkIfExists: true)
    } else {
        ref_ch = 'random' //how to send a non-path null along for auto-selection of file downstream?
    }
    inp_ch = Channel.fromPath(params.inpath, checkIfExists: true)
    out_ch = Channel.fromPath(params.outpath, checkIfExists: true)

    FIND_INFILES(
        inp_ch
    )

    INFILE_HANDLING(
        FIND_INFILES.out.found_infiles,
        inp_ch,
        out_ch,
        ref_ch
    )

    RUN_PARSNP(
        INFILE_HANDLING.out.handled_infiles,
        INFILE_HANDLING.out.ref_path,
        INFILE_HANDLING.out.tmp_path,
        out_ch
    )

    EXTRACT_SNPS(
        RUN_PARSNP.out.ran_parsnp,
        out_ch
    )

    PAIRWISE_DISTANCES(
        EXTRACT_SNPS.out.extracted_snps,
        out_ch
    )

    DISTANCE_MATRIX(
        PAIRWISE_DISTANCES.out.calculated_snp_distances,
        out_ch
    )

    if (params.recombination != false) {
        EXTRACT_FASTA(
            RUN_PARSNP.out.parsnp_xmfa,
            out_ch
        )
    }

    if (params.recombination == 'gubbins') {
        INFER_RECOMBINATION_GUBBINS(
            EXTRACT_FASTA.out.parsnp_fasta,
            RUN_PARSNP.out.parsnp_tree,
            out_ch
        )
    } else if (params.recombination == 'cfml') {
        INFER_RECOMBINATION_CFML(
            EXTRACT_FASTA.out.parsnp_fasta,
            RUN_PARSNP.out.parsnp_tree,
            out_ch
        )
    } else if (params.recombination == 'both') {
        INFER_RECOMBINATION_GUBBINS(
            EXTRACT_FASTA.out.parsnp_fasta,
            RUN_PARSNP.out.parsnp_tree,
            out_ch
        )
        INFER_RECOMBINATION_CFML(
            EXTRACT_FASTA.out.parsnp_fasta,
            RUN_PARSNP.out.parsnp_tree,
            out_ch
        )
    }

//     MASK_RECOMBINATION(
//         EXTRACT_FASTA.out.parsnp_fasta,
//         RUN_PARSNP.out.parsnp_tree,
//         INFER_RECOMBINATION_GUBBINS.out.recombination_positions,
//         'gubbins'
//         out_ch
//     )
//
//     REINFER_TREE(
//         MASK_RECOMBINATION.out.masked_fasta,
//         out_ch
//     )
}


/*
==============================================================================
                        Completion e-mail and summary                         
==============================================================================
*/
workflow.onComplete {
    workDir = new File("${workflow.workDir}")

    println """
    Pipeline Execution Summary
    --------------------------
    Workflow Version : ${workflow.version}
    Nextflow Version : ${nextflow.version}
    Command Line     : ${workflow.commandLine}
    Resumed          : ${workflow.resume}
    Completed At     : ${workflow.complete}
    Duration         : ${workflow.duration}
    Success          : ${workflow.success}
    Exit Code        : ${workflow.exitStatus}
    Error Report     : ${workflow.errorReport ?: '-'}
    Launch Dir       : ${workflow.launchDir}
    """
}

workflow.onError {
    def err_msg = """\
        Error summary
        ---------------------------
        Completed at: ${workflow.complete}
        exit status : ${workflow.exitStatus}
        workDir     : ${workflow.workDir}
        
        ??? extra error messages to include ???
        """
        .stripIndent()
/*    sendMail(
        to: "${USER}@cdc.gov",
        subject: 'workflow error',
        body: err_msg,
        charset: UTF-8
        // attach: '/path/stderr.log.txt'
    )
*/
}
