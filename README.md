# Assembly SNPs Workflow


## Steps in the workflow
1. Have conda or docker installed
2. Identify all FastA assembly files in a given input path
    - recognized extensions are:  fa, fasta, fas, fna, fsa, fa.gz, fasta.gz, fas.gz, fna.gz, fsa.gz
3. Create a tmp dir and symlink to the assembly files with file extensions removed from all sample names
4. Run parsnp, which generates a phylogenetic tree "parsnp.tree"
    - optionally specify the reference file (`--reference <FILE>`); otherwise largest filesize is chosen automatically as the reference file
5. Extract a FastA file "SNPs.fa" with harvesttools of all samples with only the SNP positions
6. Perform all pairwise comparisons to tabulate the number of SNPs each sample pair have between them
7. Generate a matrix table summarizing pairwise distances between all samples
8. Optionally, classify SNPs as being due to recombination

![workflow](images/workflow_v1.0.0.png)


## Example install
```
# install
git clone git@github.com:chrisgulvik/wf-assembly-snps.git $HOME/wf-assembly-snps
cd $HOME/wf-assembly-snps
```

## Run with conda
```
# make conda and nextflow available for use
module load conda nextflow
# run workflow
nextflow run main.nf --outpath OUTPATH_DIR --inpath INPUT_DIR -with-dag flow.png
```

## Run with docker
```
# make sure docker engine is running, then
# run workflow
nextflow run -profile docker main.nf --outpath OUTPATH_DIR --inpath INPUT_DIR -with-dag flow.png
# run workflow with recombination (options are 'gubbins', 'cfml', or 'both')
nextflow run -profile docker main.nf --outpath OUTPATH_DIR --inpath INPUT_DIR --recombination gubbins -with-dag flow.png
```

## Workflow output
```
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
|-- .log
|   |-- stderr.nextflow.txt
|   `-- stdout.nextflow.txt
|-- clonalframeml
|   |-- clonalframeml.importation_status.txt
|   |-- clonalframeml.labelled_tree.newick
|   |-- clonalframeml_masked_recombination.fasta
|   `-- clonalframeml_masked_recombination.tree
|-- gubbins
|   |-- gubbins.node_labelled.final_tree.tre
|   |-- gubbins.recombination_predictions.gff
|   |-- gubbins_masked_recombination.fasta
|   `-- gubbins_masked_recombination.tree
|-- parsnp
|   |-- SNP-distances.matrix.tsv
|   |-- SNP-distances.pairs.tsv
|   |-- SNPs.fa
|   |-- parsnp.fasta
|   |-- parsnp.ggr
|   |-- parsnp.tree
|   `-- parsnp.xmfa
|-- pipeline_dag.2022-02-18\ 12:52:15.svg
|-- report.2022-02-18\ 12:52:15.html
|-- timeline.2022-02-18\ 12:52:15.html
`-- trace.2022-02-18\ 12:52:15.txt

# cleanup
rm -rf .nextflow .nextflow.log* work/ OUTPATH_DIR/
```


## Misc Notes
- `-resume` uses cached results
- `-stub` uses example output files for any process with an uncommented stub block
- `-with-dag dag.png` to make dag

- scripts needs to be in ./bin for nextflow to be able to find them

- doesn't seem possible to tell nextflow where to find conda, it only checks your path

- how to stop appending -ue to bash
    - add `process.shell = ['/bin/bash']` to nextflow.config

