process MASK_RECOMBINATION {

    publishDir "${params.outpath}/${recombination_method}",
        mode: "${params.publish_dir_mode}",
        pattern: "*.fasta"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}" }

    container "snads/mask-recombination@sha256:0df4f5e26b2beeb2a180c2e4d75148cde55d4cc62585b5053d6606c6623d33e4"

    input:
        path alignment
        path node_labelled_tree
        path recombination_positions
        val  recombination_method

    output:
        path("${recombination_method}_masked_recombination.fasta"), emit: masked_alignment
        path(".command.out")
        path(".command.err")
        path "versions.yml",                                        emit: versions

    shell:
    '''
    mask_recombination.py \
        --alignment !{alignment} \
        --format !{recombination_method} \
        --rec_positions !{recombination_positions} \
        --tree !{node_labelled_tree} \
        > !{recombination_method}_masked_recombination.fasta

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        python: $(python --version)
    END_VERSIONS
    '''

}