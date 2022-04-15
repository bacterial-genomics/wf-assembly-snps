process PARSNP {

    publishDir "${params.outpath}",
        mode: "${params.publish_dir_mode}",
        pattern: "parsnp/*"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}" }

    label "process_medium"
    label "error_retry"

    // TODO: replace conda and null below with appropriate conda env, singularity image for this process
    //conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        null :
        'snads/parsnp@sha256:0dc552de1cf91fb30aa25eb09b4d9eef841abae989760c937b0006dacd165377' }"

    input:
    path reference
    path tmpdir

    output:
    path "parsnp/parsnp.ggr",  emit: ggr
    path "parsnp/parsnp.xmfa", emit: xmfa
    path "parsnp/parsnp.tree", emit: tree
    path ".command.out"
    path ".command.err"
    path "versions.yml",       emit: versions

    script:
    """
    source bash_functions.sh

    # TEMPORARY hack to avoid the v1.2 issue that thinks input ref file has
    #  aligned nucleotide sequences if the deflines contain hyphens
    sed -i "s/-//g" ${reference}

    # Perform the parsnp system call
    if [[ "${params.curated_input}" == "true" && "${params.tree_method}" == "fasttree" ]]; then
        msg "INFO: starting parsnp system call with curated input (will retain all sequences)"
        parsnp -v -d ${tmpdir} -r ${reference} -o parsnp -p ${task.cpus} -P ${params.max_partition_size} --verbose -c --use-fasttree
    elif [[ "${params.curated_input}" == "false" && "${params.tree_method}" == "fasttree" ]]; then
        msg "INFO: starting parsnp system call with uncurated input (will check sequence similarity)"
        parsnp -v -d ${tmpdir} -r ${reference} -o parsnp -p ${task.cpus} -P ${params.max_partition_size} --verbose --use-fasttree
    elif [[ "${params.curated_input}" == "true" && "${params.tree_method}" == "raxml" ]]; then
        msg "INFO: starting parsnp system call with curated input (will retain all sequences)"
        parsnp -v -d ${tmpdir} -r ${reference} -o parsnp -p ${task.cpus} -P ${params.max_partition_size} --verbose -c
    elif [[ "${params.curated_input}" == "false" && "${params.tree_method}" == "raxml" ]]; then
        msg "INFO: starting parsnp system call with uncurated input (will check sequence similarity)"
        parsnp -v -d ${tmpdir} -r ${reference} -o parsnp -p ${task.cpus} -P ${params.max_partition_size} --verbose
    fi

    msg "INFO: finished parsnp system call"

    # Verify output
    if ! verify_file_minimum_size "parsnp/parsnp.ggr" 'gingr' '${params.min_ggr_size}'; then
      echo "ERROR: Parsnp failed, ggr file is too small." >&2
      exit 1
    elif ! verify_file_minimum_size "parsnp/parsnp.xmfa" 'alignment' '${params.min_xmfa_size}'; then
      echo "ERROR: Parsnp failed, alignment file is too small. Did you use a curated input directory with any too-unrelated genomes?" >&2
      exit 1
    fi

    # Remove reference label from tip in tree file
    sed -i "s/.ref//1" parsnp/parsnp.tree

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parsnp: \$(parsnp --version | sed 's/^/    /')
    END_VERSIONS
    """

}
