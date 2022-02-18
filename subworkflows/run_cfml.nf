include {
    INFER_RECOMBINATION_CFML;
    MASK_RECOMBINATION;
    REINFER_TREE
} from '../modules/assembly-snps.nf'


workflow RUN_CFML {
    recombination_method = Channel.from('clonalframeml')
    take:
        parsnp_fasta
        parsnp_tree
    main:
        INFER_RECOMBINATION_CFML(
            parsnp_fasta,
            parsnp_tree
        )
        MASK_RECOMBINATION(
            parsnp_fasta,
            INFER_RECOMBINATION_CFML.out.node_labelled_tree,
            INFER_RECOMBINATION_CFML.out.recombination_positions,
            recombination_method
        )
        REINFER_TREE(
            MASK_RECOMBINATION.out.masked_fasta,
            recombination_method
        )
}