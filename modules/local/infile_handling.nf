process INFILE_HANDLING {

    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}" }

    // TODO: replace conda and null below with appropriate conda env, singularity image for this process
    //conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "ubuntu:focal"

    input:
    path input
    path reference

    output:
    path "ref/*",        emit: reference
    path "tmp",          emit: tmpdir
    path "tmp/*"
    path ".command.out"
    path ".command.err"
    path "versions.yml", emit: versions

    shell:
    '''
    source bash_functions.sh

    shopt -s nullglob
    compressed_asm=( "!{input}"/*.{fa,fas,fsa,fna,fasta}.gz )
    plaintext_asm=( "!{input}"/*.{fa,fas,fsa,fna,fasta} )
    msg "INFO: ${#compressed_asm[@]} compressed assemblies found"
    msg "INFO: ${#plaintext_asm[@]} plain text assemblies found"
    shopt -u nullglob

    # Create a new input path only containing assembly files, with each
    #  filename lacking an extension. This avoids file extensions from showing
    #  up in the sample names in parsnp output files.
    mkdir ./tmp
    for file in "${compressed_asm[@]}" "${plaintext_asm[@]}"; do
      if [ "${file: -3}" == '.gz' ]; then
        base="$(basename "${file}" .gz)"
        pref_no_ext="${base%.*}"
        # Store decompressed assemblies in outpath area to avoid write
        #  permission issues in the inpath. Also remove file extensions.
        gunzip -c "${file}" > ./tmp/"${pref_no_ext}"
      else
        base="$(basename "${file}")"
        pref_no_ext="${base%.*}"
        ln -sf "$(pwd)"/"${file}" ./tmp/"${pref_no_ext}"
      fi
    done

    # Confirm >1 file exists in the newly formed tmp/ subdirectory
    shopt -s nullglob
    tmp_asm=( ./tmp/* )
    shopt -u nullglob
    if [[ ${#tmp_asm[@]} -le 1 ]]; then
      msg "ERROR: unable to create tmp subdir to store input assemblies" >&2
      exit 1
    fi
    msg "INFO: created tmp path with all ${#tmp_asm[@]} input assembly files"

    # Parsnp is reference-based. Handle a specified reference file, otherwise
    #  use the largest filesize for reproducibility (not randomly selected).
    mkdir ./ref
    if [ "!{reference}" = 'largest' ]; then
      # Get paths to real input files, following symlinks
      ls -l tmp/* | awk '{print $NF}' > to_sort.txt
      # Get name of largest file
      largest_size=$(du -a `cat to_sort.txt` | sort -n | tail -n 1 | awk '{print $2}' | xargs basename)
      # Remove extension if largest file was symlinked (to match name in tmp/)
      largest_size="$(echo $largest_size | sed -E -e "s/\\.fa|\\.fas|\\.fsa|\\.fna|\\.fasta$//")"
      reference=./ref/"${largest_size}"
      if [[ -L "./tmp/"${largest_size}"" ]]; then
        # Copy reference if it is a symlink to input dir
        cp "$(readlink ./tmp/"${largest_size}")" "${reference}"
        rm "./tmp/"${largest_size}""
      else
        # Move reference if it is a real file in tmp/
        mv -f ./tmp/"${largest_size}" "${reference}"
      fi
      msg "INFO: ${largest_size} autoselected as reference with largest filesize"
    else
      # Handle if the reference file is compressed and in a different input path
      #  than the specified --inpath <PATH> which would have already been handled.
      if [[ !{reference} == *.gz ]]; then
        base="$(basename "!{reference}" .gz)"
        pref_no_ext="${base%.*}"
        rm -f ./tmp/"${pref_no_ext}" # handles whether it's there or not
        gunzip -c "!{reference}" > ./ref/"${pref_no_ext}"
        msg "INFO: specified compressed reference !{reference} copied to ref"
      else
        base="$(basename "!{reference}")"
        pref_no_ext="${base%.*}"
        rm -f ./tmp/"${pref_no_ext}" # handles whether it's there or not
        cp -L "!{reference}" ./ref/"${pref_no_ext}"
        msg "INFO: specified plaintext reference !{reference} copied to ref"
      fi
    fi

    # Replace problematic characters in file names
    shopt -s nullglob  # removes literal pattern as match if no real files match
    for f in ./tmp/*"|"*; do
        mv "$f" "${f//|/_}"
    done
    for f in ./ref/*"|"*; do
        mv "$f" "${f//|/_}"
    done
    shopt -u nullglob

    # Confirm exactly 1 reference file exists
    shopt -s nullglob
    ref_asm=( ./ref/* )
    shopt -u nullglob
    if [[ ${#ref_asm[@]} -ne 1 ]]; then
      msg "ERROR: unable to create the ref/<FILE> reference file" >&2
      exit 1
    fi
    msg "INFO: will use the ${ref_asm[0]} reference file"

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        ubuntu: $(cat /etc/issue)
    END_VERSIONS
    '''

}
