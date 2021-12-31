#!/usr/bin/env python3

# This script finds pairwise distances of entries in an alignment file
# TODO: add option to read in from stdin
# TODO: add estimate
# TODO: how to handle ambig characters?

import sys
import argparse
from Bio import SeqIO
from multiprocessing import Pool
# from simulate import *
import functools
import numpy as np

def main():
    args = parse_args()

    if args.simulate:
        if args.verbose: print("Simulating an alignment.")
        args.infile = "./sim_alignment.fasta"
        simulate_alignment(sim_outfile = args.infile)

    # Read in sequences
    if args.verbose: print("Reading in sequences.")
    with open(args.infile) as f:
        n_seqs = int(len(f.readlines()) / 2)
    seqs = [ '' for i in range(n_seqs) ]
    i = 0
    for rec in SeqIO.parse(args.infile, "fasta"):
        seqs[i] = rec.seq
        i += 1

    # Work through tuples of sequence idxs, calculating distances between respective sequences
    dist_matrix = [ [ 0 for i in range(n_seqs) ] for j in range(n_seqs) ]
    jobs = [ (i, j) for i in range(n_seqs) for j in range(i + 1, n_seqs) ]
    if args.verbose: print("Calculating " + str(len(jobs)) + " pairwise distances.")
    with Pool(args.numcpus) as p:
        results = p.map(functools.partial(pairwise_distance, seqs), jobs)

    # Unpack results list into upper-triangular distance matrix, write out
    if args.verbose: print("Reformatting into upper-triangular distance matrix.")
    for k in range(len(jobs)):
        dist_matrix[jobs[k][0]][jobs[k][1]] = results[k]

    if args.verbose: print("Writing out distance matrix.")
    np.savetxt(args.outfile, dist_matrix, fmt='%i')


def parse_args():
    parser = argparse.ArgumentParser(description = "Find pairwise distances of entries in an alignment file.")
    parser.add_argument('infile')
    parser.add_argument('outfile', nargs="?", type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('-n', '--numcpus', type=int, default=1, help="Number of CPUs (default: 1).")
    parser.add_argument('-e', '--estimate', action='store_true', help="Estimate the number of pairwise distances using random sampling. 1/4 of all pairwise bases will be analyzed instead of 100%%.")
    parser.add_argument('-s', '--estfreq', type=float, default=0.25, help="(to be used with -e) The frequency at which to analyze positions for pairwise differences.")
    parser.add_argument('-t', '--simulate', action='store_true', help="Simulate an alignment, disregard infile. For testing only.")
    parser.add_argument('-v', '--verbose', action='store_true', help="Print progress updates.")
    return parser.parse_args()


def pairwise_distance(seqs, idx_tuple):
    seq_i = seqs[idx_tuple[0]]
    seq_j = seqs[idx_tuple[1]]
    dist = 0
    for nt_pair in zip(seq_i, seq_j):
        if nt_pair[0] != nt_pair[1]:
            dist += 1
    return dist


if __name__ == '__main__':
    main()
