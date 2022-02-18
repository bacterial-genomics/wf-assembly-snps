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

![workflow](images/workflow_v1.0.0.svg)

## Requirements
* Nextflow
* Docker or Singularity

## Install
```
git clone git@github.com:chrisgulvik/wf-assembly-snps.git
cd wf-assembly-snps
```

[//]: # (## Run with conda)

[//]: # (```)

[//]: # (# make conda and nextflow available for use)

[//]: # (module load conda nextflow)

[//]: # (# run workflow)

[//]: # (nextflow run main.nf --outpath OUTPATH_DIR --inpath INPUT_DIR -with-dag flow.png)

[//]: # (```)

## Set-up for Aspen cluster
``` 
# Add these Singularity variables to $HOME/.bashrc
SINGULARITY_BASE=/scicomp/scratch/$USER
export SINGULARITY_TMPDIR=$SINGULARITY_BASE/singularity.tmp
export SINGULARITY_CACHEDIR=$SINGULARITY_BASE/singularity.cache
mkdir -pv $SINGULARITY_TMPDIR $SINGULARITY_CACHEDIR

# Restart session

module load nextflow/21.04.3

# Run workflow with -profile singularity
```

## Run workflow
```
# With Docker
nextflow run -profile docker main.nf --outpath OUTPATH_DIR --inpath INPUT_DIR
# With Singularity
nextflow run -profile singularity main.nf --outpath OUTPATH_DIR --inpath INPUT_DIR
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


## Misc User Notes
- `-resume` uses cached results
- `-stub` uses example output files for any process with an uncommented stub block

## Misc Dev Notes
- scripts needs to be in ./bin for nextflow to be able to find them
- doesn't seem possible to tell nextflow where to find conda, it only checks your path
- how to stop appending -ue to bash
    - add `process.shell = ['/bin/bash']` to nextflow.config

