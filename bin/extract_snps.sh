#!/usr/bin/env bash


parsnp_ggr=$1
snps_outfile=$2

source bash_functions.sh

# Extract only SNP positions from the gingr file
msg "INFO: starting SNP extraction"
harvesttools \
  -i "$parsnp_ggr" \
  -S "$snps_outfile"
msg "INFO: finished SNP extraction"

if ! verify_file_minimum_size "${snps_outfile}" 'SNP' '10c'; then
  msg "ERROR: harvesttools failed to extract SNPs" >&2
  exit 1
fi
sed -i 's/\.ref//1' "${snps_outfile}"
