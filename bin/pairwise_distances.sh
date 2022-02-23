#!/usr/bin/env bash

cpus=$1
snps_file=$2

source bash_functions.sh

# Calculate pairwise SNP distances from all sample pairs
snp_distances="SNP-distances.pairs.tsv"
msg "INFO: starting to calculate pairwise SNP distances"
#pairwiseDistances.pl \
# -n "${cpus}" \
# "${snps_file}" \
# | sort -k3,3n \
# > "${snp_distances}"
pairwiseDistances.py \
  -n "${cpus}" \
  "${snps_file}" \
   | sort -k3,3n \
   > "${snp_distances}"
msg "INFO: finished calculating pairwise SNP distances"

if ! verify_file_minimum_size "${snp_distances}" 'SNP distances' '20c'; then
  msg "ERROR: pairwise SNP distances not calculated" >&2
  exit 1
fi
