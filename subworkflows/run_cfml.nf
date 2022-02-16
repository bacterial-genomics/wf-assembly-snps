include {
    INFER_RECOMBINATION_CFML;
    MASK_RECOMBINATION;
    REINFER_TREE
} from '../modules/assembly-snps.nf'

workflow RUN_CFML {
    take:
        parsnp_fasta
        parsnp_tree
        out_ch
    main:
        INFER_RECOMBINATION_CFML(
            parsnp_fasta,
            parsnp_tree,
            out_ch
        )
        MASK_RECOMBINATION(
            parsnp_fasta,
            parsnp_tree,
            INFER_RECOMBINATION_CFML.out.recombination_positions,
            channel.from('clonalframeml'),
            out_ch
        )
        REINFER_TREE(
            MASK_RECOMBINATION.out.masked_fasta,
            channel.from('clonalframeml'),
            out_ch
        )
}