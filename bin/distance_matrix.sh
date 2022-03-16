#!/usr/bin/env bash


snp_distances=$1

source bash_functions.sh

# Create a 2-dimensional matrix from pairwise SNP distances of all sample pairs
snp_matrix=SNP-distances.matrix.tsv

msg "INFO: starting to form SNP matrix"
pairwiseTo2d.py \
 -i "${snp_distances}" \
 -o "${snp_matrix}" \
 --sort
msg "INFO: finished SNP matrix formation"

if ! verify_file_minimum_size "${snp_matrix}" 'SNP matrix' '20c'; then
  msg "ERROR: SNP distance matrix not formed" >&2
  exit 1
fi
sed -i "s/\t-/\t0/g" "${snp_matrix}"
