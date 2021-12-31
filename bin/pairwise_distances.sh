#!/usr/bin/env bash


output_path=$1
cpus=$2

source bash_functions.sh

# Calculate pairwise SNP distances from all sample pairs
snps_outfile="${output_path}"/SNPs.fa
snp_distances="${output_path}"/SNP-distances.pairs.tsv
msg "INFO: starting to calculate pairwise SNP distances"
pairwiseDistances.pl \
 -n "${cpus}" \
 "${snps_outfile}" \
 | sort -k3,3n \
 > "${snp_distances}"
# calc_hamming_dist.py \
#  -n "${cpus}" \
#  "${snp_distances}" \
#  "${snps_outfile}"
msg "INFO: finished calculating pairwise SNP distances"

if ! verify_file_minimum_size "${snp_distances}" 'SNP distances' '20c'; then
  msg "ERROR: pairwise SNP distances not calculated" >&2
  exit 1
fi

touch "pairwise_distances.success.txt"
