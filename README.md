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
# install
git clone git@github.com:chrisgulvik/wf-assembly-snps.git $HOME/wf-assembly-snps
cd $HOME/wf-assembly-snps

# make conda and nextflow available for use
module load conda nextflow

# run the workflow with defaults
nextflow run main.nf --outpath OUTPATH_DIR --inpath INPUT_DIR -with-dag flow.png

# view final output file
cat OUTPATH_DIR/SNP-distances.matrix.tsv
-   16-090  16-100  16-127  16-146  16-151  16-155
16-090  0   31  24  7   32  35
16-100  31  0   33  32  3   6
16-127  24  33  0   25  34  37
16-146  7   32  25  0   33  36
16-151  32  3   34  33  0   3
16-155  35  6   37  36  3   0

# view final output dir structure
tree -a OUTPATH_DIR/
OUTPATH_DIR/
├── .log
│   ├── stderr.nextflow.txt
│   └── stdout.nextflow.txt
├── parsnpAligner.log.gz
├── parsnp.ggr
├── parsnp.tree
├── SNP-distances.matrix.tsv
├── SNP-distances.pairs.tsv
├── SNPs.fa.gz
└── trace.2021-36-30 22:36:26.txt

# cleanup
rm -rf .nextflow .nextflow.log* work/ OUTPATH_DIR/
```


## Misc Notes
`-resume` uses cached results
`-with-dag dag.png` to make dag

scripts needs to be in ./bin for nextflow to be able to find them

doesn't seem possible to tell nextflow where to find conda,
it only checks your path

- how to stop appending -ue to bash
    - add `process.shell = ['/bin/bash']` to nextflow.config
