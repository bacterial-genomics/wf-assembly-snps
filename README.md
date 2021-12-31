# Assembly SNPs Workflow


## Steps in the workflow
1. Identify all FastA assembly files in a given input path
    - recognized extensions are:  fa, fasta, fas, fna, fsa, fa.gz, fasta.gz, fas.gz, fna.gz, fsa.gz
2. Create a tmp dir and symlink to the assembly files with file extensions removed from all sample names
3. Run parsnp, which generates a phylogenetic tree "parsnp.tree"
    - optionally specify the reference file (`--reference <FILE>`); otherwise largest filesize is chosen automatically as the reference file
4. Extract a FastA file "SNPs.fa" with harvesttools of all samples with only the SNP positions
5. Perform all pairwise comparisons to tabulate the number of SNPs each sample pair have between them
6. Generate a matrix table summarizing pairwise distances between all samples
7. Optionally, classify SNPs as being due to recombination

![workflow](images/workflow_v1.0.0.png)


## Example install and run
```
git clone git@github.com:chrisgulvik/wf-assembly-snps.git $HOME
cd $HOME/wf-assembly-snps
nextflow run main.nf --outpath OUTPATH_DIR --inpath INPUT_DIR -with-dag flow.png
```


## Misc Notes
`-resume` uses cached results
`-with-dag dag.png` to make dag

scripts needs to be in ./bin for nextflow to be able to find them

doesn't seem possible to tell nextflow where to find conda,
it only checks your path

- how to stop appending -ue to bash
    - add `process.shell = ['/bin/bash']` to nextflow.config
