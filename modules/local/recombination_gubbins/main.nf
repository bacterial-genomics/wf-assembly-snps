process RECOMBINATION_GUBBINS {

    tag { "${meta.snp_package}" }
    label "process_medium"
    container "snads/gubbins@sha256:391a980312096f96d976f4be668d4dea7dda13115db004a50e49762accc0ec62"

    input:
    tuple val(meta), path(core_alignment_fasta)
    tuple val(meta_alignment), path(alignment_files)

    output:
    tuple val(meta), path("*.{txt,tree}"), emit: positions_and_tree
    path(".command.{out,err}")
    path("versions.yml")                 , emit: versions

    shell:
    '''
    source bash_functions.sh

    msg "INFO: Performing recombination using Gubbins."

    run_gubbins.py \
      --starting-tree "!{meta.snp_package}.tree" \
      --prefix "!{meta.snp_package}-Gubbins" \
      "!{core_alignment_fasta}"

    # Rename output files
    mv "!{meta.snp_package}-Gubbins.recombination_predictions.gff" \
      "!{meta.snp_package}-Gubbins.recombination_positions.txt"

    mv "!{meta.snp_package}-Gubbins.node_labelled.final_tree.tre" \
      "!{meta.snp_package}-Gubbins.labelled_tree.tree"

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        gubbins: $(run_gubbins.py --version | sed 's/^/    /')
    END_VERSIONS
    '''
}
