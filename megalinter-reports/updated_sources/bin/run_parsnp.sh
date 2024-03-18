#!/usr/bin/env bash

input_path=$1
output_path=$2
reference_path=$3
cpus=$4

source bash_functions.sh

command -v parsnp >/dev/null 2>&1 || {
  echo 'ERROR: parsnp not found' >&2
  exit 1
}

# TEMPORARY hack to avoid the v1.2 issue that thinks input ref file has
#  aligned nucleotide sequences if the deflines contain hyphens
sed -i 's/-//g' "${reference_path}"

# Log binary path and the version used
msg "Parsnp path $(command -v parsnp)"
msg "Parsnp version $(parsnp --version 2>/dev/null | sed 's/parsnp //1')"

# Perform the parsnp system call
msg "INFO: starting parsnp system call"
parsnp \
  -v \
  -c \
  -d "${input_path}" \
  -r "${reference_path}" \
  -o "${output_path}" \
  -p "${cpus}"
msg "INFO: finished parsnp system call"

if ! verify_file_minimum_size "${output_path}/parsnp.ggr" 'gingr' '100k'; then
  echo """ERROR: Parsnp failed.
  Check ${output_path}/parsnpAligner.log""" >&2
  exit 1
fi
rm -rf "${output_path}"/.{ref,tmp}
sed -i 's/\.ref//1' "${output_path}/parsnp.tree"

touch "run_parsnp.success.txt"
