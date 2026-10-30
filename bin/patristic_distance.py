#!/usr/bin/env python3

import argparse
import csv
from itertools import combinations

from Bio import Phylo


def get_leaf_names(tree):
    """Return terminal node names from a tree, preserving tree order."""
    leaves = [leaf.name for leaf in tree.get_terminals()]

    if any(name is None for name in leaves):
        raise ValueError("All terminal nodes in the tree must have names.")

    if len(set(leaves)) != len(leaves):
        raise ValueError("Terminal node names must be unique.")

    return leaves


def write_all_vs_all_distances(tree, leaf_names, outfile):
    """Write unique pairwise patristic distances as a sparse TSV."""
    with open(outfile, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["genome1", "genome2", "patristic_dist"])

        for genome1, genome2 in combinations(leaf_names, 2):
            distance = float(tree.distance(genome1, genome2))
            writer.writerow([genome1, genome2, format(distance, ".7f")])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="python3 patristic_distance.py",
        description=(
            "Generate an all-vs-all sparse TSV of patristic distances from a "
            "NEWICK tree."
        ),
    )

    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Path to NEWICK file",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="./patristic_distances.tsv",
        help="Path and name of output TSV",
    )
    args = parser.parse_args()

    tree = Phylo.read(args.input, "newick")
    leaf_names = get_leaf_names(tree)
    write_all_vs_all_distances(tree, leaf_names, args.output)
