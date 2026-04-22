/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Defaults
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_assemblysnps_pipeline'

// Added
include { QUAST } from '../modules/nf-core/quast/main'
include { GUBBINS } from '../modules/nf-core/gubbins/main'
include { PARSNP } from '../modules/local/parsnp/main'
include { FASTTREE } from '../modules/nf-core/fasttree/main'
include { SNPSITES } from '../modules/nf-core/snpsites/main'
include { SNPDISTS } from '../modules/nf-core/snpdists/main'
include { CLONALFRAMEML } from '../modules/nf-core/clonalframeml/main'
include { IQTREE } from '../modules/nf-core/iqtree/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ASSEMBLYSNPS {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    //
    // PREPROCESSING
    //

    // Filter assemblies by size
    ch_samplesheet
        .branch { meta, fasta ->
            pass: fasta[0].size() > params.min_fasta_size
            fail: fasta[0].size() <= params.min_fasta_size
        }
        .set { ch_samplesheet_filtered }

    ch_samplesheet_filtered.fail.tap { ch_samplesheet_fail_log }
    ch_samplesheet_fail_log.view( it -> "SAMPLE FAIL: Length of ${it[0].id} < ${params.min_fasta_size} bytes" )

    ch_reference = file(params.reference, checkIfExists: true)

    //
    // MODULE: QUAST
    //

    QUAST (
        ch_samplesheet_filtered.pass,
        [[],[]], // tuple val(meta2), path(fasta)
        [[],[]], // tuple val(meta3), path(gff)
    )
    ch_quast_multiqc = QUAST.out.results

    //
    // MODULE: Parsnp
    //

    ch_parsnp = ch_samplesheet_filtered.pass
        .map { meta, fasta -> fasta }
        .collect()

    PARSNP (
        ch_parsnp,
        ch_reference
    )

    //
    // RECOMBINATION DETECTION
    //

    if (params.run_gubbins) {

        //
        // MODULE: GUBBINS
        //

        ch_gubbins = PARSNP.out.aln
        ch_tree = IQTREE.out.phylogeny.map { meta, tree -> tree }

        GUBBINS ( ch_gubbins, ch_tree )
    }

    if (params.run_clonalframeml) {

        //
        // MODULE: ClonalFrameML
        //

        ch_clonalframeml = FASTTREE.out.phylogeny
            .combine( SNPSITES.out.fasta )
            .map { newick, msa -> [[], newick, msa] }

        CLONALFRAMEML (
            ch_clonalframeml
        )
    }

    //
    // MODULE: SnpSites
    //

    ch_snpsites = PARSNP.out.aln

    SNPSITES (
        ch_snpsites
    )

    //
    // MODULE: SNPdists
    //

    ch_snpdists = SNPSITES.out.fasta.map { aln -> [[], aln] }

    SNPDISTS (
        ch_snpdists
    )

    //
    // MODULE: FastTree
    //

    ch_fasttree = SNPSITES.out.fasta

    FASTTREE (
        ch_fasttree
    )

    //
    // MODULE: IQ-TREE
    //

    ch_iqtree = SNPSITES.out.fasta
        .map { aln -> [[], aln, []] }

    IQTREE (
        ch_iqtree,
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        []
    )

    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'assemblysnps_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
