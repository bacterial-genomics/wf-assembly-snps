include {
    INFER_RECOMBINATION_GUBBINS;
    MASK_RECOMBINATION;
    REINFER_TREE
} from '../modules/assembly-snps.nf'


workflow RUN_GUBBINS {
    recombination_method = Channel.from('gubbins')
    take:
        parsnp_fasta
        parsnp_tree
    main:
        INFER_RECOMBINATION_GUBBINS(
            parsnp_fasta,
            parsnp_tree
        )
        MASK_RECOMBINATION(
            parsnp_fasta,
            INFER_RECOMBINATION_GUBBINS.out.node_labelled_tree,
            INFER_RECOMBINATION_GUBBINS.out.recombination_positions,
            recombination_method
        )
        REINFER_TREE(
            MASK_RECOMBINATION.out.masked_fasta,
            recombination_method
        )
}