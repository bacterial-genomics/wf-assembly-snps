process CORE_GENOME_ALIGNMENT_PARSNP {

    tag { "${meta_input.snp_package}" }
    label "process_medium"
    container "quay.io/biocontainers/parsnp@sha256:b46999fb9842f183443dd6226b461c1d8074d4c1391c1f2b1e51cc20cee8f8b2"

    input:
    tuple val(meta_input)    , path("genomes/")
    tuple val(meta_reference), path(reference_file)

    output:
    tuple val(meta_input), path("Parsnp_Output_File.tsv"), emit: qc_filecheck
    tuple val(meta_input), path("*{ggr,xmfa,tree,fa.gz}"), emit: output
    path(".command.{out,err}")
    path("versions.yml")                                 , emit: versions

    shell:
    curatedInput = params.curated_input               ? "--curated"      : ""
    treeMethod   = (params.tree_method == "fasttree") ? "--use-fasttree" : ""
    '''
    source bash_functions.sh

    # TEMPORARY hack to avoid the v1.2 issue that thinks input ref file has
    #  aligned nucleotide sequences if the deflines contain hyphens
    sed -i "s/-//g" !{reference_file}

    msg "INFO: Performing the parsnp system call."

    parsnp \
      --sequences genomes/ \
      --reference !{reference_file} \
      --output-dir Parsnp \
      --verbose \
      !{treeMethod} \
      !{curatedInput} \
      --threads !{task.cpus} \
      --max-partition-size !{params.max_partition_size}

    # Move files to current directory
    mv Parsnp/parsnp.ggr Parsnp.ggr
    mv Parsnp/parsnp.xmfa Parsnp.xmfa
    mv Parsnp/parsnp.snps.mblocks Parsnp.SNPs.fa

    if [[ !{params.tree_method} == "fasttree" ]]; then
      mv Parsnp/log/fasttree.out Parsnp.tree
    else
      mv Parsnp/parsnp.tree Parsnp.tree
    fi

    # Verify output
    echo -e "Sample name\tQC step\tOutcome (Pass/Fail)" > Parsnp_Output_File.tsv
    for file in Parsnp.ggr Parsnp.xmfa Parsnp.SNPs.fa; do
      if verify_minimum_file_size "${file}" "${file} Output File" "!{params.min_parsnp_output_filesize}"; then
        echo -e "NaN\t${file} Output File\tPASS" >> Parsnp_Output_File.tsv
      else
        echo -e "NaN\t${file} Output File\tFAIL" >> Parsnp_Output_File.tsv
      fi
    done

    # Remove the 1 additional suffix Parsnp adds to the reference sample `.ref`
    if [[ $(grep -o -n '.ref' Parsnp.tree | wc -l) -eq 1 ]]; then
      sed -i 's/.ref//1' Parsnp.tree Parsnp.SNPs.fa
    fi

    gzip -9f Parsnp.SNPs.fa

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        parsnp: $(parsnp --version | sed 's/^/    /')
    END_VERSIONS
    '''
}
