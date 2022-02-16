#!/usr/bin/env bash


input_path=$1
output_path=$2
reference_path=$3

source bash_functions.sh

shopt -s nullglob
compressed_asm=( "${input_path}"/*.{fa,fas,fsa,fna,fasta}.gz )
plaintext_asm=( "${input_path}"/*.{fa,fas,fsa,fna,fasta} )
shopt -u nullglob

# Create a new input path only containing assembly files, with each
#  filename lacking an extension. This avoids file extensions from showing
#  up in the sample names in parsnp output files.
mkdir "${output_path}"/.tmp
for file in "${compressed_asm[@]}" "${plaintext_asm[@]}"; do
  if [ "${file: -3}" == '.gz' ]; then
    base="$(basename "${file}" .gz)"
    pref_no_ext="${base%.*}"
    # Store decompressed assemblies in outpath area to avoid write
    #  permission issues in the inpath. Also remove file extensions.
    gunzip -c "${file}" > "${output_path}"/.tmp/"${pref_no_ext}"
  else
    base="$(basename "${file}")"
    pref_no_ext="${base%.*}"  
    # Save space by just making symlinks with filenames lacking file
    #  extensions.
    ln -sf "${file}" "${output_path}"/.tmp/"${pref_no_ext}"
  fi
done

# Confirm >1 file exists in the newly formed .tmp/ subdirectory
shopt -s nullglob
tmp_asm=( "${output_path}"/.tmp/* )
shopt -u nullglob
if [[ ${#tmp_asm[@]} -le 1 ]]; then
  msg "ERROR: unable to create .tmp subdir to store input assemblies" >&2
  exit 1
fi
msg "INFO: created .tmp path with all ${#tmp_asm[@]} input assembly files"

# Parsnp is reference-based. Handle a specified reference file, otherwise
#  use the largest filesize for reproducibility (not randomly selected).
mkdir "${output_path}"/.ref
if [ "${reference_path}" == 'largest' ]; then
  largest_size=$(ls -SLr "${output_path}"/.tmp/* | tail -n 1 | xargs basename)
  reference_path="${output_path}"/.ref/"${largest_size}"
  mv -f "${output_path}"/.tmp/"${largest_size}" "${reference_path}"
  msg "INFO: ${largest_size} autoselected as reference with largest filesize"
else
  # Handle if the reference file is compressed and in a different input path
  #  than the specified --inpath <PATH> which would have already been handled.
  if [ "${reference_path: -3}" == '.gz' ]; then
    base="$(basename "${reference_path}" .gz)"
    pref_no_ext="${base%.*}"
    rm -f "${output_path}"/.tmp/"${pref_no_ext}" # handles whether it's there or not
    gunzip -c "${reference_path}" > "${output_path}"/.ref/"${pref_no_ext}"
    msg "INFO: specified compressed reference ${reference_path} copied to .ref"
  else
    base="$(basename "${reference_path}")"
    pref_no_ext="${base%.*}"
    rm -f "${output_path}"/.tmp/"${pref_no_ext}" # handles whether it's there or not
    cp -L "${reference_path}" "${output_path}"/.ref/"${pref_no_ext}"
    msg "INFO: specified plaintext reference ${reference_path} copied to .ref"
  fi
fi

# Confirm exactly 1 reference file exists
shopt -s nullglob
ref_asm=( "${output_path}"/.ref/* )
shopt -u nullglob
if [[ ${#ref_asm[@]} -ne 1 ]]; then
  msg "ERROR: unable to create the .ref/<FILE> reference file" >&2
  exit 1
fi
msg "INFO: will use the ${ref_asm[0]} reference file"

touch "infile_handling.success.txt"
