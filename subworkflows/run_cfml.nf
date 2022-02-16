include {
    INFER_RECOMBINATION_CFML;
    MASK_RECOMBINATION;
    REINFER_TREE
} from '../modules/assembly-snps.nf'


workflow RUN_CFML {
    take:
        extract_fasta_success
        parsnp_fasta
        parsnp_tree
        outpath
    main:
        INFER_RECOMBINATION_CFML(
            parsnp_fasta,
            parsnp_tree,
            outpath
        )
        MASK_RECOMBINATION(
            INFER_RECOMBINATION_CFML.out.infer_recombination_success,
            parsnp_fasta,
            INFER_RECOMBINATION_CFML.out.node_labelled_tree,
            INFER_RECOMBINATION_CFML.out.recombination_positions,
            channel.from('clonalframeml'),
            outpath
        )
        REINFER_TREE(
            MASK_RECOMBINATION.out.mask_recombination_success,
            MASK_RECOMBINATION.out.masked_fasta,
            channel.from('clonalframeml'),
            outpath
        )
}