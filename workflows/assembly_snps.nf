/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowSNPS.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet or directory not specified!' }
if (params.ref) { ch_ref_input = file(params.ref) } else { ch_ref_input = [] }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// CONFIGS: Import configs for this workflow
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES: Local modules
//
include { INFILE_HANDLING_UNIX                             } from "../modules/local/infile_handling_unix/main"
include { INFILE_HANDLING_UNIX as REF_INFILE_HANDLING_UNIX } from "../modules/local/infile_handling_unix/main"

include { CORE_GENOME_ALIGNMENT_PARSNP                     } from "../modules/local/core_genome_alignment_parsnp/main"
include { CONVERT_GINGR_TO_FASTA_HARVESTTOOLS              } from "../modules/local/convert_gingr_to_fasta_harvesttools/main"

include { CALCULATE_PAIRWISE_DISTANCES_SNP_DISTS           } from "../modules/local/calculate_pairwise_distances_snp_dists/main"
include { CREATE_SNP_DISTANCE_MATRIX_SNP_DISTS             } from "../modules/local/create_snp_distance_matrix_snp_dists/main"
include { MASK_RECOMBINANT_POSITIONS_BIOPYTHON             } from "../modules/local/mask_recombinant_positions_biopython/main"
include { CREATE_MASKED_SNP_DISTANCE_MATRIX_SNP_DISTS      } from "../modules/local/create_masked_snp_distance_matrix_snp_dists/main"

include { BUILD_PHYLOGENETIC_TREE_PARSNP                   } from "../modules/local/build_phylogenetic_tree_parsnp/main"

include { CONVERT_TSV_TO_EXCEL_PYTHON                      } from "../modules/local/convert_tsv_to_excel_python/main"
include { CREATE_EXCEL_RUN_SUMMARY_PYTHON                  } from "../modules/local/create_excel_run_summary_python/main"

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK                                      } from "../subworkflows/local/input_check"
include { INPUT_CHECK as REF_INPUT_CHECK                   } from "../subworkflows/local/input_check"
include { RECOMBINATION                                    } from "../subworkflows/local/recombination"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CREATE CHANNELS FOR INPUT PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if ( toLower(params.snp_package) == "parsnp" ) {
    ch_snp_package = "Parsnp"
} else {
    ch_snp_package = "Parsnp"
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Convert input to lowercase
def toLower(it) {
    it.toString().toLowerCase()
}

// Check QC filechecks for a failure
def qcfilecheck(process, qcfile, inputfile) {
    qcfile.map{ meta, file -> [ meta, [file] ] }
            .join(inputfile)
            .map{ meta, qc, input ->
                data = []
                qc.flatten().each{ data += it.readLines() }

                if ( data.any{ it.contains('FAIL') } ) {
                    line = data.last().split('\t')
                    if (line.first() != "NaN") {
                        log.warn("${line[1]} QC check failed during process ${process} for sample ${line.first()}")
                    } else {
                        log.warn("${line[1]} QC check failed during process ${process}")
                    }
                } else {
                    [ meta, input ]
                }
            }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ASSEMBLY_SNPS {

    // SETUP: Define empty channels to concatenate certain outputs
    ch_versions             = Channel.empty()
    ch_qc_filecheck         = Channel.empty()
    ch_output_summary_files = Channel.empty()

    /*
    ================================================================================
                            Preprocess input data
    ================================================================================
    */

    // SUBWORKFLOW: Check input for samplesheet or pull inputs from directory
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // Check input files meet size criteria
    INFILE_HANDLING_UNIX (
        INPUT_CHECK.out.input_files
    )
    ch_versions     = ch_versions.mix(INFILE_HANDLING_UNIX.out.versions)
    ch_qc_filecheck = ch_qc_filecheck.concat(INFILE_HANDLING_UNIX.out.qc_filecheck)

    // qcfilecheck function is not needed for this channel due to implementation within module
    ch_input_files  = INFILE_HANDLING_UNIX.out.input_files
                        .collect()
                        .map{
                            def meta = [:]
                            meta['snp_package'] = ch_snp_package
                            [ meta, it ]
                        }

    // SUBWORKFLOW: Check reference input for samplesheet or pull inputs from directory
    if (params.ref) {
        REF_INPUT_CHECK (
            ch_ref_input
        )
        ch_versions = ch_versions.mix(REF_INPUT_CHECK.out.versions)

        // Check input files meet size criteria
        REF_INFILE_HANDLING_UNIX (
            REF_INPUT_CHECK.out.input_files
        )
        ch_versions        = ch_versions.mix(REF_INFILE_HANDLING_UNIX.out.versions)
        ch_qc_filecheck    = ch_qc_filecheck.concat(REF_INFILE_HANDLING_UNIX.out.qc_filecheck)

        // qcfilecheck function is not needed for this channel due to implementation within module
        ch_reference_files = REF_INFILE_HANDLING_UNIX.out.input_files
                                .collect()
                                .map{
                                    def meta = [:]
                                    meta['snp_package'] = ch_snp_package
                                    [ meta, it ]
                                }

    } else {
        // Reverse sort items by size and use first item as reference and remove it from ch_input_files
        ch_reference_files = ch_input_files
                                .map{
                                    meta, file ->
                                        def file_new = file.collect().sort{ it.size() }.reverse()
                                        [ meta, file_new[0] ]
                                }.collect()

        ch_input_files     = ch_input_files
                                .map{
                                    meta, file ->
                                    def file_new = file.collect().sort{ it.size() }.reverse()
                                    file_new.remove(file_new[0])
                                    [ meta, file_new ]
                                }.collect()
    }

    // Check for reference file in input channel; if so, remove it
    ch_input_files.join(ch_reference_files)
        .map{
            meta, input, ref ->
                if ( input.contains(ref) ) {
                    input.remove(ref)
                    [ meta, input ]
                } else {
                    [ meta, input ]
                }
        }
        .set { ch_input_files }

    /*
    ================================================================================
                            Perform core alignment
    ================================================================================
    */

    if ( toLower(params.snp_package) == "parsnp" ) {
        // PROCESS: Run ParSNP to generate core genome alignment and phylogeny
        CORE_GENOME_ALIGNMENT_PARSNP (
            ch_input_files,
            ch_reference_files
        )
        ch_versions        = ch_versions.mix(CORE_GENOME_ALIGNMENT_PARSNP.out.versions)
        ch_qc_filecheck    = ch_qc_filecheck.concat(CORE_GENOME_ALIGNMENT_PARSNP.out.qc_filecheck)

        ch_alignment_files = qcfilecheck(
                                "CORE_GENOME_ALIGNMENT_PARSNP",
                                CORE_GENOME_ALIGNMENT_PARSNP.out.qc_filecheck,
                                CORE_GENOME_ALIGNMENT_PARSNP.out.output
                            )

        // PROCESS: Convert Parsnp Gingr output file to FastA format for recombination
        CONVERT_GINGR_TO_FASTA_HARVESTTOOLS (
            ch_alignment_files
        )
        ch_versions = ch_versions.mix(CONVERT_GINGR_TO_FASTA_HARVESTTOOLS.out.versions)
        ch_qc_filecheck = ch_qc_filecheck.concat(CONVERT_GINGR_TO_FASTA_HARVESTTOOLS.out.qc_filecheck)

        ch_core_alignment_fasta = qcfilecheck(
                                    "CONVERT_GINGR_TO_FASTA_HARVESTTOOLS",
                                    CONVERT_GINGR_TO_FASTA_HARVESTTOOLS.out.qc_filecheck,
                                    CONVERT_GINGR_TO_FASTA_HARVESTTOOLS.out.core_alignment
                                )
    }

    /*
    ================================================================================
                            Calculate distances
    ================================================================================
    */

    // PROCESS: Calculate pairwise genome distances
    CALCULATE_PAIRWISE_DISTANCES_SNP_DISTS (
        ch_alignment_files
    )
    ch_versions             = ch_versions.mix(CALCULATE_PAIRWISE_DISTANCES_SNP_DISTS.out.versions)
    ch_output_summary_files = ch_output_summary_files.mix(CALCULATE_PAIRWISE_DISTANCES_SNP_DISTS.out.snp_distances.map{ meta, file -> file })

    // PROCESS: Reformat pairwise genome distances into matrix
    CREATE_SNP_DISTANCE_MATRIX_SNP_DISTS (
        ch_alignment_files
    )
    ch_versions             = ch_versions.mix(CREATE_SNP_DISTANCE_MATRIX_SNP_DISTS.out.versions)
    ch_output_summary_files = ch_output_summary_files.mix(CREATE_SNP_DISTANCE_MATRIX_SNP_DISTS.out.distance_matrix.map{ meta, file -> file })

    /*
    ================================================================================
                            Identify and mask recombinant positions
    ================================================================================
    */

    // SUBWORKFLOW: Infer SNPs due to recombination
    RECOMBINATION (
        ch_core_alignment_fasta,
        ch_alignment_files
    )
    ch_versions = ch_versions.mix(RECOMBINATION.out.versions)

    // PROCESS: Mask recombinant positions
    MASK_RECOMBINANT_POSITIONS_BIOPYTHON (
        RECOMBINATION.out.recombinants,
        ch_core_alignment_fasta.collect()
    )
    ch_versions = ch_versions.mix(MASK_RECOMBINANT_POSITIONS_BIOPYTHON.out.versions)

    // PROCESS: Create SNP distance matrix on masked FastA file
    CREATE_MASKED_SNP_DISTANCE_MATRIX_SNP_DISTS (
        MASK_RECOMBINANT_POSITIONS_BIOPYTHON.out.masked_alignment
    )
    ch_versions             = ch_versions.mix(CREATE_MASKED_SNP_DISTANCE_MATRIX_SNP_DISTS.out.versions)
    ch_output_summary_files = ch_output_summary_files.mix(CREATE_MASKED_SNP_DISTANCE_MATRIX_SNP_DISTS.out.distance_matrix.map{ meta, file -> file })

    /*
    ================================================================================
                            Build phylogenetic tree
    ================================================================================
    */

    // PROCESS: Infer phylogenetic tree based on masked positions
    BUILD_PHYLOGENETIC_TREE_PARSNP (
        MASK_RECOMBINANT_POSITIONS_BIOPYTHON.out.masked_alignment
    )
    ch_versions     = ch_versions.mix(BUILD_PHYLOGENETIC_TREE_PARSNP.out.versions)
    ch_qc_filecheck = ch_qc_filecheck.concat(BUILD_PHYLOGENETIC_TREE_PARSNP.out.qc_filecheck)

    ch_final_tree   = qcfilecheck(
                            "BUILD_PHYLOGENETIC_TREE_PARSNP",
                            BUILD_PHYLOGENETIC_TREE_PARSNP.out.qc_filecheck,
                            BUILD_PHYLOGENETIC_TREE_PARSNP.out.tree
                        )

    /*
    ================================================================================
                        Collect QC information
    ================================================================================
    */

    // Collect QC file check information
    ch_qc_filecheck = ch_qc_filecheck
                        .map{ meta, file -> file }
                        .collectFile(
                            name:       "Summary.QC_File_Checks.tsv",
                            keepHeader: true,
                            storeDir:   "${params.outdir}/Summaries",
                            sort:       'index'
                        )

    ch_output_summary_files = ch_output_summary_files.mix(ch_qc_filecheck.collect())

    /*
    ================================================================================
                        Convert TSV outputs to Excel XLSX
    ================================================================================
    */

    if (params.create_excel_outputs) {
        CREATE_EXCEL_RUN_SUMMARY_PYTHON (
            ch_output_summary_files.collect()
        )
        ch_versions = ch_versions.mix(CREATE_EXCEL_RUN_SUMMARY_PYTHON.out.versions)

        CONVERT_TSV_TO_EXCEL_PYTHON (
            CREATE_EXCEL_RUN_SUMMARY_PYTHON.out.summary
        )
        ch_versions = ch_versions.mix(CONVERT_TSV_TO_EXCEL_PYTHON.out.versions)
    }

    /*
    ================================================================================
                        Collect version information
    ================================================================================
    */

    // Collect version information
    ch_versions
        .unique()
        .collectFile(
            name:     "software_versions.yml",
            storeDir: params.tracedir
        )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
