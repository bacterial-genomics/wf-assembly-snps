#!/usr/bin/env python3

'''
Mask recombination positions identified using Gubbins or ClonalFrameML.

Original author: Sarah Nadeau
Link to original script: https://github.com/SarahNadeau/docker-builds/blob/sarah-dev/mask_recombination/1.0/mask_recombination.py
'''

import argparse
from Bio import SeqIO, Phylo
from Bio.Seq import MutableSeq
from functools import partial
import gzip
import re
import sys


def parseArgs():
    parser = argparse.ArgumentParser(
        description="Mask SNPs due to recombination in an alignment file."
    )
    parser.add_argument(
        "--alignment",
        help="FastA-formatted core genome alignment. Optionally gzipped with suffix '.gz'.",
        required=True
    )
    parser.add_argument(
        "--rec_positions",
        help="ClonalFrameML or Gubbins output: positions affected by recombination.",
        required=True
    )
    parser.add_argument(
        "--tree",
        help="ClonalFrameML or Gubbins output: phylogeny with internal nodes labelled.",
        required=True
    )
    parser.add_argument(
        "--format",
        help="Program used to estimate recombination positions and phylogeny: 'clonalframeml' or 'gubbins'.",
        default='clonalframeml'
    )
    parser.add_argument(
        "--verbose",
        help="Log intermediate results.",
        action='store_true',
        default=False
    )
    return parser.parse_args()


def main():
    opt = parseArgs()

    logmsg("Reading in phylogeny, recombination position inputs.")
    tree = Phylo.read(opt.tree, 'newick')
    if opt.format == 'clonalframeml':
        recombinations = get_recombinations_clonalframeml(opt.rec_positions)
    elif opt.format == 'gubbins':
        recombinations = get_recombinations_gubbins(opt.rec_positions)
    else:
        raise ValueError("--format must be either 'clonalframeml' or 'gubbins'")
    if opt.verbose:
        logmsg("Recombinations:")
        logmsg(recombinations)
    validate_recombinations(recombinations, tree)

    logmsg("Traversing tree to get all tips to apply mask to for each recombination event.")
    masks = get_masks(tree, recombinations)
    if opt.verbose:
        logmsg("Masks:")
        logmsg(masks)

    logmsg("Iterating through sequence file, masking recombination positions.")
    mask_alignment(opt.alignment, masks)


def logmsg(*opt, **kwopt):
    print(*opt, file=sys.stderr, **kwopt)


# Parse ClonalFrameML recombination estimates
# Returns dictionary with node name keys and recombination position values
def get_recombinations_clonalframeml(file):
    recombinations = {}
    f = open(file)
    for line in f:
        if not re.match("^Node.*Beg.*End$", line):
            parsed_line = line.strip().split()
            node_name = parsed_line[0]
            start = int(parsed_line[1]) - 1  # ClonalFrameML results are 1-indexed
            stop = int(parsed_line[2])  # ClonalFrameML results are end-inclusive
            if node_name in recombinations.keys():
                recombinations[node_name] = recombinations[node_name] + [(start, stop)]
            else:
                recombinations[node_name] = [(start, stop)]
    f.close()
    return recombinations


# Parse Gubbins recombination estimates
# Returns dictionary with node name keys and recombination position values
def get_recombinations_gubbins(file):
    recombinations = {}
    f = open(file)
    for line in f:
        if not re.match("^##", line):
            parsed_line = line.strip().split()
            node_name = parsed_line[8].split("\";")[0].split("->")[1]
            start = int(parsed_line[3]) - 1  # Gubbins results are 1-indexed
            stop = int(parsed_line[4])  # Gubbins results are end-inclusive
            if node_name in recombinations.keys():
                recombinations[node_name] = recombinations[node_name] + [(start, stop)]
            else:
                recombinations[node_name] = [(start, stop)]
    f.close()
    return recombinations


# Check that all nodes in recombination file are present in the phylogeny
def validate_recombinations(recombinations, tree):
    nonterminal_nodes = [node.name for node in tree.get_nonterminals()]
    terminal_nodes = [node.name for node in tree.get_terminals()]
    tree_nodes = set(nonterminal_nodes + terminal_nodes)
    recombination_nodes = {str(key) for key, value in recombinations.items()}
    missing_tree_nodes = recombination_nodes - tree_nodes
    if len(missing_tree_nodes) > 0:
        raise AssertionError(
            "Some nodes in the recombination file do not appear in the phylogeny: {}".format(missing_tree_nodes))


# Summarize recombination positions at tips in tree
# Returns dictionary with tip name keys and recombination position values
def get_masks(tree, recombinations):
    masks = {}
    for node_name in recombinations.keys():
        matches = tree.find_clades(name=node_name)
        first_match = next(matches)
        tip_names = [tip.name for tip in first_match.get_terminals()]
        for tip_name in tip_names:
            if tip_name in masks.keys():
                masks[tip_name] = masks[tip_name] + recombinations[node_name]
            else:
                masks[tip_name] = recombinations[node_name]
    return masks


# Writes out a FastA file with recombination positions masked with 'N'
def mask_alignment(alignment_file, masks):
    _open = partial(gzip.open, mode='rt') if alignment_file.endswith('.gz') else open
    with _open(alignment_file) as f:
        first_seq = True
        for rec in SeqIO.parse(f, "fasta"):
            seq_len = len(rec.seq)
            if first_seq:
                aln_len = seq_len
                first_seq = False
            if len(rec.seq) != aln_len:  # Input must be aligned
                raise AssertionError(
                    "Some sequences appear un-aligned! " +
                    "{} is length {} instead of expected {}.".format(rec.id, len(rec.seq), aln_len))
            if rec.id in masks.keys():
                new_seq = MutableSeq(rec.seq)
                for mask in masks[rec.id]:
                    if max(mask) > seq_len:
                        raise ValueError("A recombination event spans the end of the alignment.")
                    start = mask[0]
                    end = mask[1]
                    new_seq[start:end] = 'N' * (end - start)
            else:
                new_seq = rec.seq
            print(">" + rec.id)
            print(new_seq)


if __name__ == "__main__":
    main()
