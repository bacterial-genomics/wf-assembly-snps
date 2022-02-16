include {
    INFER_RECOMBINATION_GUBBINS;
    MASK_RECOMBINATION;
    REINFER_TREE
} from '../modules/assembly-snps.nf'


workflow RUN_GUBBINS {
    take:
        extract_fasta_success
        parsnp_fasta
        parsnp_tree
        outpath
    main:
        INFER_RECOMBINATION_GUBBINS(
            parsnp_fasta,
            parsnp_tree,
            outpath
        )
        MASK_RECOMBINATION(
            INFER_RECOMBINATION_GUBBINS.out.infer_recombination_success,
            parsnp_fasta,
            parsnp_tree,
            INFER_RECOMBINATION_GUBBINS.out.recombination_positions,
            channel.from('gubbins'),
            outpath
        )
        REINFER_TREE(
            MASK_RECOMBINATION.out.mask_recombination_success,
            MASK_RECOMBINATION.out.masked_fasta,
            channel.from('gubbins'),
            outpath
        )
}