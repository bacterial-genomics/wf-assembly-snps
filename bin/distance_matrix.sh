#!/usr/bin/env bash


output_path=$1

source bash_functions.sh

# Create a 2-dimensional matrix from pairwise SNP distances of all sample pairs
snp_distances="${output_path}"/SNP-distances.pairs.tsv
snp_matrix="${output_path}"/SNP-distances.matrix.tsv

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

# Cleanup
rm -f "${output_path}"/{all_mumi.ini,parsnpAligner.ini,parsnp.{rec,xmfa},psnn.ini}
rm -f "${output_path}"/*.ref
rm -rf "${output_path}"/{tmp,.tmp,.ref}
pigz -9f "${output_path}"/{SNPs.fa,parsnpAligner.log}

touch "distance_matrix.success.txt"
