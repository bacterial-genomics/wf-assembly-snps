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
include { EXTRACT_SNP_POSITIONS_PARSNP                     } from "../modules/local/extract_snp_positions_parsnp/main"
include { CALCULATE_PAIRWISE_DISTANCES_BIOPYTHON           } from "../modules/local/calculate_pairwise_distances_biopython/main"
include { CREATE_SNP_DISTANCE_MATRIX_BIOPYTHON             } from "../modules/local/create_snp_distance_matrix_biopython/main"
include { CONVERT_XMFA_FASTA_PYTHON                        } from "../modules/local/convert_xmfa_fasta_python/main"
include { MASK_RECOMBINANT_POSITIONS_BIOPYTHON             } from "../modules/local/mask_recombinant_positions_biopython/main"

include { BUILD_PHYLOGENETIC_TREE_PARSNP                   } from "../modules/local/build_phylogenetic_tree_parsnp/main"

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

if (params.ref) {
    ch_ref_input = Channel.fromPath(params.ref, checkIfExists: true)
} else {
    ch_ref_input = Channel.empty()
}

if ( toLower(params.aligner) == "parsnp" ) {
    ch_aligner = "Parsnp"
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Convert params.aai to lowercase
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
                            meta['aligner'] = ch_aligner
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
                                    meta['aligner'] = ch_aligner
                                    [ meta, it ]
                                }

    } else {
        // Use first item as reference and remove it from ch_input_files
        ch_reference_files = ch_input_files
                                .map{
                                    meta, file ->
                                        def file_new = file.collect().sort()
                                        [ meta, file_new[0] ]
                                }.collect()

        ch_input_files     = ch_input_files
                                .map{
                                    meta, file ->
                                    def file_new = file.collect().sort()
                                    file_new.remove(file_new[0])
                                    [ meta, file_new ]
                                }.collect()
    }

    /*
    ================================================================================
                            Perform core alignment
    ================================================================================
    */

    if ( toLower(params.aligner) == "parsnp" ) {
        // PROCESS: Run ParSNP to generate core genome alignment and phylogeny
        CORE_GENOME_ALIGNMENT_PARSNP (
            ch_input_files,
            ch_reference_files
        )
        ch_versions        = ch_versions.mix(CORE_GENOME_ALIGNMENT_PARSNP.out.versions)
        ch_qc_filecheck    = ch_qc_filecheck.concat(CORE_GENOME_ALIGNMENT_PARSNP.out.qc_filecheck)

        ch_gingr_alignment = qcfilecheck(
                                "CORE_GENOME_ALIGNMENT_PARSNP",
                                CORE_GENOME_ALIGNMENT_PARSNP.out.qc_filecheck,
                                CORE_GENOME_ALIGNMENT_PARSNP.out.gingr_alignment
                            )

        ch_core_alignment  = qcfilecheck(
                                "CORE_GENOME_ALIGNMENT_PARSNP",
                                CORE_GENOME_ALIGNMENT_PARSNP.out.qc_filecheck,
                                CORE_GENOME_ALIGNMENT_PARSNP.out.core_alignment
                            )

        ch_parsnp_snps     = qcfilecheck(
                                "CORE_GENOME_ALIGNMENT_PARSNP",
                                CORE_GENOME_ALIGNMENT_PARSNP.out.qc_filecheck,
                                CORE_GENOME_ALIGNMENT_PARSNP.out.snps
                            )
    }
    /*
    ================================================================================
                            Calculate distances
    ================================================================================
    */

    // PROCESS: Get SNP alignment from ParSNP .ggr file
    EXTRACT_SNP_POSITIONS_PARSNP (
        ch_gingr_alignment
    )
    ch_versions       = ch_versions.mix(EXTRACT_SNP_POSITIONS_PARSNP.out.versions)
    ch_qc_filecheck   = ch_qc_filecheck.concat(EXTRACT_SNP_POSITIONS_PARSNP.out.qc_filecheck)

    ch_extracted_snps = qcfilecheck(
                            "EXTRACT_SNP_POSITIONS_PARSNP",
                            EXTRACT_SNP_POSITIONS_PARSNP.out.qc_filecheck,
                            EXTRACT_SNP_POSITIONS_PARSNP.out.snps
                        )

    // PROCESS: Calculate pairwise genome distances
    CALCULATE_PAIRWISE_DISTANCES_BIOPYTHON (
        ch_extracted_snps
    )
    ch_versions      = ch_versions.mix(CALCULATE_PAIRWISE_DISTANCES_BIOPYTHON.out.versions)
    ch_qc_filecheck  = ch_qc_filecheck.concat(CALCULATE_PAIRWISE_DISTANCES_BIOPYTHON.out.qc_filecheck)

    ch_snp_distances = qcfilecheck(
                            "CALCULATE_PAIRWISE_DISTANCES_BIOPYTHON",
                            CALCULATE_PAIRWISE_DISTANCES_BIOPYTHON.out.qc_filecheck,
                            CALCULATE_PAIRWISE_DISTANCES_BIOPYTHON.out.snp_distances
                        )

    // PROCESS: Reformat pairwise genome distances into matrix
    CREATE_SNP_DISTANCE_MATRIX_BIOPYTHON (
        CALCULATE_PAIRWISE_DISTANCES_BIOPYTHON.out.snp_distances
    )
    ch_versions            = ch_versions.mix(CREATE_SNP_DISTANCE_MATRIX_BIOPYTHON.out.versions)
    ch_qc_filecheck        = ch_qc_filecheck.concat(CREATE_SNP_DISTANCE_MATRIX_BIOPYTHON.out.qc_filecheck)

    ch_snp_distance_matrix = qcfilecheck(
                                "CREATE_SNP_DISTANCE_MATRIX_BIOPYTHON",
                                CREATE_SNP_DISTANCE_MATRIX_BIOPYTHON.out.qc_filecheck,
                                CREATE_SNP_DISTANCE_MATRIX_BIOPYTHON.out.distance_matrix
                            )

    /*
    ================================================================================
                            Identify and mask recombinant positions
    ================================================================================
    */

    // SUBWORKFLOW: Infer SNPs due to recombination, mask them, re-infer phylogeny
    RECOMBINATION (
        ch_parsnp_snps,
        CORE_GENOME_ALIGNMENT_PARSNP.out.phylogeny
    )
    ch_versions = ch_versions.mix(RECOMBINATION.out.versions)

    // PROCESS: Mask recombinant positions
    MASK_RECOMBINANT_POSITIONS_BIOPYTHON (
        RECOMBINATION.out.recombinants,
        ch_parsnp_snps.collect()
    )
    ch_versions = ch_versions.mix(MASK_RECOMBINANT_POSITIONS_BIOPYTHON.out.versions)

    /*
    ================================================================================
                            Build phylogenetic tree
    ================================================================================
    */

    // PROCESS: Build phylogenetic tree
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
                        Collect version and QC information
    ================================================================================
    */

    // Collect version information
    ch_versions
        .unique()
        .collectFile(
            name:     "software_versions.yml",
            storeDir: params.tracedir
        )

    // Collect QC file check information
    ch_qc_filecheck = ch_qc_filecheck
                        .map{ meta, file -> file }
                        .collectFile(
                            name:       "Summary.QC_File_Checks.tsv",
                            keepHeader: true,
                            storeDir:   "${params.outdir}/Summaries",
                            sort:       'index'
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
