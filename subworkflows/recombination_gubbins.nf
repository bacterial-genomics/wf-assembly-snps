nextflow.enable.dsl=2

include { GUBBINS } from '../modules/local/gubbins'
include { MASK_RECOMBINATION } from '../modules/local/mask_recombination'
include { INFER_TREE } from '../modules/local/infer_tree'


workflow RECOMBINATION_GUBBINS {
    recombination_method = Channel.from('gubbins')
    take:
        alignment
        starting_tree
    main:
        GUBBINS(
            alignment,
            starting_tree
        )
        MASK_RECOMBINATION(
            alignment,
            GUBBINS.out.node_labelled_tree,
            GUBBINS.out.recombination_positions,
            recombination_method
        )
        INFER_TREE(
            MASK_RECOMBINATION.out.masked_alignment,
            recombination_method,
            params.tree_method
        )
}