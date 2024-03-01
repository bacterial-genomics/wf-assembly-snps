process INFILE_HANDLING_UNIX {

    tag( "${meta.id}" )
    container "ubuntu:jammy"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("${meta.id}.Initial_Input_File.tsv"), emit: qc_filecheck
    path("inputfiles/*")                                      , emit: input_files
    path(".command.{out,err}")
    path("versions.yml")                                      , emit: versions

    shell:
    // Rename files with meta.id (has spaces and periods removed)
    gzip_compressed = input.toString().contains('.gz') ? '.gz' : ''
    file_extension  = input.toString().split('.gz')[0].split('\\.')[-1]
    '''
    source bash_functions.sh

    # Rename input files to prefix and move to inputfiles dir
    mkdir inputfiles
    cp !{input} inputfiles/"!{meta.id}.!{file_extension}!{gzip_compressed}"

    # gunzip all files that end in .{gz,Gz,GZ,gZ}
    find -L inputfiles/ -type f -name '*.[gG][zZ]' -exec gunzip -f {} +

    # Filter out small inputfiles
    msg "Checking input file sizes.."
    echo -e "Sample name\tQC step\tOutcome (Pass/Fail)" > "!{meta.id}.Initial_Input_File.tsv"
    for file in inputfiles/*; do
      if verify_minimum_file_size "${file}" 'Input' "!{params.min_input_filesize}"; then
        echo -e "$(basename ${file%%.*})\tInput File\tPASS" \
        >> "!{meta.id}.Initial_Input_File.tsv"
      else
        echo -e "$(basename ${file%%.*})\tInput File\tFAIL" \
        >> "!{meta.id}.Initial_Input_File.tsv"

        rm ${file}
      fi
    done

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
      ubuntu: $(awk -F ' ' '{print $2,$3}' /etc/issue | tr -d '\\n')
    END_VERSIONS
    '''
}
