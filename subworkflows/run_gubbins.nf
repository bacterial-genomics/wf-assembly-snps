include {
    INFER_RECOMBINATION_GUBBINS;
    MASK_RECOMBINATION;
    REINFER_TREE
} from '../modules/assembly-snps.nf'

workflow RUN_GUBBINS {
    take:
        parsnp_fasta
        parsnp_tree
        out_ch
    main:
        INFER_RECOMBINATION_GUBBINS(
            parsnp_fasta,
            parsnp_tree,
            out_ch
        )
        MASK_RECOMBINATION(
            parsnp_fasta,
            parsnp_tree,
            INFER_RECOMBINATION_GUBBINS.out.recombination_positions,
            channel.from('gubbins'),
            out_ch
        )
        REINFER_TREE(
            MASK_RECOMBINATION.out.masked_fasta,
            channel.from('gubbins'),
            out_ch
        )
}