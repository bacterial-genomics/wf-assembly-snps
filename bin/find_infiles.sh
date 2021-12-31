#!/usr/bin/env bash


input_path=$1

source bash_functions.sh

shopt -s nullglob
compressed_asm=( $input_path/*.{fa,fas,fsa,fna,fasta}.gz )
plaintext_asm=( $input_path/*.{fa,fas,fsa,fna,fasta} )
shopt -u nullglob

msg "INFO: ${#compressed_asm[@]} compressed assemblies found"
msg "INFO: ${#plaintext_asm[@]} plain text assemblies found"

combo=$(( ${#compressed_asm[@]} + ${#plaintext_asm[@]} ))
if [[ ${combo} -eq 0 ]]; then
  msg "ERROR: unable to find assembly files in ${input_path}" >&2
  exit 1
fi

touch "find_infiles.success.txt"
