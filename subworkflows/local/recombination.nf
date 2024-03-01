//
// Identify recombinant positions
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES: Local modules
//
include { RECOMBINATION_GUBBINS       } from "../../modules/local/recombination_gubbins/main"
include { RECOMBINATION_CLONALFRAMEML } from "../../modules/local/recombination_clonalframeml/main"

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Convert params.aai to lowercase
def toLower(it) {
    it.toString().toLowerCase()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN RECOMBINATION WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RECOMBINATION {

    take:
    core_alignment
    phylogeny

    main:
    ch_versions      = Channel.empty()
    ch_recombination = Channel.empty()

    // Perform recombination
    if ( toLower(params.recombination) == "gubbins" ) {
        // PROCESS: Perform recombination using Gubbins
        RECOMBINATION_GUBBINS (
            core_alignment,
            phylogeny
        )
        ch_versions      = ch_versions.mix(RECOMBINATION_GUBBINS.out.versions)
        ch_recombination = RECOMBINATION_GUBBINS.out.positions_and_tree

    } else if ( toLower(params.recombination) == "clonalframeml" ) {
        // PROCESS: Perform recombination using ClonalFrameML
        RECOMBINATION_CLONALFRAMEML (
            core_alignment,
            phylogeny
        )
        ch_versions      = ch_versions.mix(RECOMBINATION_CLONALFRAMEML.out.versions)
        ch_recombination = RECOMBINATION_CLONALFRAMEML.out.positions_and_tree

    } else if ( toLower(params.recombination) == "both" ) {
        // PROCESS: Perform recombination using Gubbins
        RECOMBINATION_GUBBINS (
            core_alignment,
            phylogeny
        )
        ch_versions      = ch_versions.mix(RECOMBINATION_GUBBINS.out.versions)
        ch_recombination = ch_recombination.mix(RECOMBINATION_GUBBINS.out.positions_and_tree)

        // PROCESS: Perform recombination using ClonalFrameML
        RECOMBINATION_CLONALFRAMEML (
            core_alignment,
            phylogeny
        )
        ch_versions      = ch_versions.mix(RECOMBINATION_CLONALFRAMEML.out.versions)
        ch_recombination = ch_recombination.mix(RECOMBINATION_CLONALFRAMEML.out.positions_and_tree)
    }

    emit:
    versions     = ch_versions
    recombinants = ch_recombination
}
