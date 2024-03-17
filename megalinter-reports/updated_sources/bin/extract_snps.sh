#!/usr/bin/env bash

output_path=$1

source bash_functions.sh

# Log binary path and the version used
msg "Harvesttools path $(command -v harvesttools)"
msg "Harvesttools version $(harvesttools --version)"

# Extract only SNP positions from the gingr file
snps_outfile="${output_path}"/SNPs.fa
msg "INFO: starting SNP extraction"
harvesttools \
  -i "${output_path}"/parsnp.ggr \
  -S "${snps_outfile}"
msg "INFO: finished SNP extraction"

if ! verify_file_minimum_size "${snps_outfile}" 'SNP' '10c'; then
  msg "ERROR: harvesttools failed to extract SNPs" >&2
  exit 1
fi
sed -i 's/\.ref//1' "${snps_outfile}"

touch "extract_snps.success.txt"
