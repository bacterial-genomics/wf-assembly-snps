nextflow.enable.dsl=2

include { CFML } from '../modules/local/cfml'
include { MASK_RECOMBINATION } from '../modules/local/mask_recombination'
include { INFER_TREE } from '../modules/local/infer_tree'


workflow RECOMBINATION_CFML {
    recombination_method = Channel.from('clonalframeml')
    take:
        alignment
        starter_tree
    main:
        CFML(
            alignment,
            starter_tree
        )
        MASK_RECOMBINATION(
            alignment,
            CFML.out.node_labelled_tree,
            CFML.out.recombination_positions,
            recombination_method
        )
        INFER_TREE(
            MASK_RECOMBINATION.out.masked_alignment,
            recombination_method,
            params.tree_method
        )
    emit:
        cfml_versions = CFML.out.versions
        mask_recombination_versions = MASK_RECOMBINATION.out.versions
        infer_tree_versions = INFER_TREE.out.versions
}