#!/usr/bin/env bash


input_path=$1
output_path=$2
reference_path=$3
cpus=$4
curated_input=$5
max_partition_size=$6

source bash_functions.sh

command -v parsnp >/dev/null 2>&1 || { echo 'ERROR: parsnp not found' >&2; exit 1; }

# TEMPORARY hack to avoid the v1.2 issue that thinks input ref file has
#  aligned nucleotide sequences if the deflines contain hyphens
sed -i 's/-//g' "${reference_path}"

# Perform the parsnp system call
if [ "${curated_input}" == "true" ]; then
  msg "INFO: starting parsnp system call with curated input (will retain all sequences)"
  parsnp \
   -v \
   -c \
   -d "${input_path}" \
   -r "${reference_path}" \
   -o "${output_path}" \
   -p "${cpus}" \
   -P "${max_partition_size}" \
   --verbose \
   --use-fasttree
else
  msg "INFO: starting parsnp system call with uncurated input (will check sequence similarity)"
  parsnp \
   -v \
   -d "${input_path}" \
   -r "${reference_path}" \
   -o "${output_path}" \
   -p "${cpus}" \
   -P "${max_partition_size}" \
   --verbose \
   --use-fasttree
fi

msg "INFO: finished parsnp system call"

if ! verify_file_minimum_size "${output_path}/parsnp.ggr" 'gingr' '100k'; then
  echo """ERROR: Parsnp failed.
  Check ${output_path}/parsnpAligner.log""" >&2
  exit 1
elif ! verify_file_minimum_size "${output_path}/parsnp.xmfa" 'alignment' '100k'; then
  echo """ERROR: Parsnp failed.
  Alignment file is too small - did you use a curated input directory with any too-unrelated genomes?
  Check ${output_path}/parsnpAligner.log""" >&2
  exit 1
fi
rm -rf "${output_path}"/.{ref,tmp}
sed -i 's/\.ref//1' "${output_path}/parsnp.tree"

touch "run_parsnp.success.txt"
