process OUTFILE_CLEANUP {

    publishDir "${params.logpath}/command_outputs",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}" }

    // TODO: replace conda and null below with appropriate conda env, singularity image for this process
    //conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        null :
        'ubuntu:focal' }"

    input:
        path(snp_distances)
        path(outpath)

    output:
        path(".command.out")
        path(".command.err")
        path "versions.yml", emit: versions

    shell:
    '''
    source bash_functions.sh

    msg "INFO: starting compressing SNPs.fa file"
    gzip -9f !{outpath}/parsnp/SNPs.fa
    #pigz -9f !{outpath}/parsnp/SNPs.fa # TODO
    msg "INFO: finished compressing SNPs.fa file"

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        ubuntu: $(cat /etc/issue)
    END_VERSIONS
    '''

}