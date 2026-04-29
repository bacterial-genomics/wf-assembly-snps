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
include { QUAST                       } from '../modules/nf-core/quast/main'
include { GUBBINS                     } from '../modules/nf-core/gubbins/main'
include { PARSNP                      } from '../modules/local/parsnp/main'
include { FASTTREE                    } from '../modules/nf-core/fasttree/main'
include { SNPSITES                    } from '../modules/nf-core/snpsites/main'
include { SNPDISTS as SNPDISTS_MATRIX } from '../modules/nf-core/snpdists/main'
include { SNPDISTS as SNPDISTS_PAIRS  } from '../modules/nf-core/snpdists/main'
include { SNPDISTS as SNPDISTS_MATRIX_RECOMB } from '../modules/nf-core/snpdists/main'
include { SNPDISTS as SNPDISTS_PAIRS_RECOMB  } from '../modules/nf-core/snpdists/main'
include { CLONALFRAMEML               } from '../modules/nf-core/clonalframeml/main'
include { IQTREE                      } from '../modules/nf-core/iqtree/main'

// local
include { GGTREE } from '../modules/local/ggtree/main'
include { CLUSTER } from '../modules/local/cluster/main'
include { GGTREE as GGTREE_RECOMB } from '../modules/local/ggtree/main'
include { CLUSTER as CLUSTER_RECOMB } from '../modules/local/cluster/main'


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
    ch_snpdists = channel.empty()
    ch_snpdists_recomb = channel.empty()
    ch_clonalframeml = channel.empty()
    ch_gubbins = channel.empty()
    ch_gubbins_tree = channel.empty()
    ch_ggtree = channel.empty()
    ch_fasttree = channel.empty()
    ch_iqtree = channel.empty()

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
    ch_samplesheet_fail_log.view { it -> "SAMPLE FAIL: Length of ${it[0].id} < ${params.min_fasta_size} bytes" }

    //ch_reference = file(params.reference, checkIfExists: true)


    //
    // MODULE: QUAST
    //

    if ( !params.skip_quast ) {
        QUAST (
            ch_samplesheet_filtered.pass,
            [[],[]], // tuple val(meta2), path(fasta)
            [[],[]], // tuple val(meta3), path(gff)
        )
        ch_quast_multiqc = QUAST.out.results
        //ch_multiqc_files = ch_multiqc_files.mix(ch_quast_multiqc)
    }

    //
    // MODULE: Parsnp
    //

    ch_parsnp = ch_samplesheet_filtered.pass
        .map { meta, fasta -> fasta }
        .collect()

    PARSNP (
        ch_parsnp,
        []
    )

    //
    // MODULE: SnpSites
    //

    ch_snpsites = PARSNP.out.aln

    SNPSITES (
        ch_snpsites
    )

    //
    // Phylogenetic tree construction
    // MODULES: FastTree, IQ-TREE
    //

    if (params.run_fasttree) {
        ch_fasttree = ch_fasttree.mix(SNPSITES.out.fasta)

        FASTTREE (
            ch_fasttree
        )

        ch_clonalframeml = FASTTREE.out.phylogeny
            .combine( SNPSITES.out.fasta )
            .map { newick, msa -> [ [], newick, msa ] }

        ch_gubbins_tree = ch_gubbins_tree.mix(FASTTREE.out.phylogeny)
        ch_ggtree = ch_ggtree
            .mix(FASTTREE.out.phylogeny)

    } else {
        ch_iqtree = ch_iqtree.mix(SNPSITES.out.fasta.map { aln -> [[ id: "iqtree"], aln, []] })

        IQTREE (
            ch_iqtree,
            [],[],[],[],[],[],[],[],[],[],[],[]
        )

        ch_gubbins_tree = ch_gubbins_tree.mix(IQTREE.out.phylogeny.map { meta, tree -> tree })

        ch_clonalframeml = ch_clonalframeml.mix(IQTREE.out.phylogeny)
            .map { meta, aln -> aln }
            .combine( SNPSITES.out.fasta )
            .map { newick, msa -> [ [], newick, msa ] }

        ch_ggtree = ch_ggtree
            .mix(IQTREE.out.phylogeny)
            .map { meta, tree -> tree }
    }

    //
    // RECOMBINATION DETECTION
    // MODULES: Gubbins, ClonalFrameML
    //

    if (params.run_gubbins) {

        ch_gubbins = ch_gubbins.mix(
            SNPSITES.out.fasta
            .collect()
            .combine ( ch_gubbins_tree.collect() )
         )

        GUBBINS ( ch_gubbins )

        ch_snpdists_recomb = ch_snpdists_recomb.mix( GUBBINS.out.fasta.map { aln -> [[], aln] } )

    } else if (params.run_clonalframeml) {
        
        CLONALFRAMEML ( ch_clonalframeml )
        
        ch_snpdists_recomb = ch_snpdists_recomb.mix( CLONALFRAMEML.out.fasta )

    }


    //
    // MODULE: SNPdists
    //

    SNPDISTS_MATRIX_RECOMB (
        ch_snpdists_recomb
    )

    SNPDISTS_PAIRS_RECOMB (
        ch_snpdists_recomb
    )

    ch_snpdists = ch_snpdists.mix( SNPSITES.out.fasta.map { aln -> [[], aln] } )

    SNPDISTS_MATRIX (
        ch_snpdists
    )

    SNPDISTS_PAIRS (
        ch_snpdists
    )

    //
    // Tree visualization and clustering
    //

    ch_cluster = SNPDISTS_MATRIX.out.tsv.map { meta, tsv -> tsv }

    CLUSTER ( ch_cluster, params.snp_threshold )

    ch_cluster_recomb = SNPDISTS_MATRIX_RECOMB.out.tsv.map { meta, tsv -> tsv }

    CLUSTER_RECOMB ( ch_cluster_recomb, params.snp_threshold )

    ch_ggtree_reg = ch_ggtree
        .combine( CLUSTER.out.clusters )

    GGTREE ( ch_ggtree_reg )

    ggtree_recomb = ch_ggtree
        .combine( CLUSTER_RECOMB.out.clusters )

    GGTREE_RECOMB ( ggtree_recomb )

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
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
