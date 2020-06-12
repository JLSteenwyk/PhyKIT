import getopt
import logging
import math
import os.path
import statistics as stat
import sys

from Bio import Phylo
import numpy as np


from .base import Tree

class InternalBranchStats(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        internal_branch_lengths, stats = self.calculate_internal_branch_stats(tree)
        if self.verbose:
            for internal_branch_length in internal_branch_lengths:
                print(internal_branch_length)
        else:
            print(f"mean: {stats['mean']}")
            print(f"median: {stats['median']}")
            print(f"25th percentile: {stats['twenty_fifth']}")
            print(f"75th percentile: {stats['seventy_fifth']}")
            print(f"minimum: {stats['minimum']}")
            print(f"maximum: {stats['maximum']}")
            print(f"standard deviation: {stats['standard_deviation']}")
            print(f"variance: {stats['variance']}")

    def process_args(self, args):
        return dict(tree_file_path=args.tree, verbose=args.verbose)

    def calculate_internal_branch_stats(self, tree):
        # save internal branch lengths to internal_branch_lengths
        internal_branch_lengths = []
        for internal_branch in tree.get_nonterminals():
            if internal_branch.branch_length != None:
                internal_branch_lengths.append(internal_branch.branch_length)

        # If the phylogeny had no branch lengths, inform user and quit
        if len(internal_branch_lengths) == 0:
            print("Calculating internal branch statistics requires a phylogeny with branch lengths.")
            sys.exit()
        
        stats = dict(
            mean               = stat.mean(internal_branch_lengths),
            median             = stat.median(internal_branch_lengths),
            twenty_fifth       = np.percentile(internal_branch_lengths, 25),
            seventy_fifth      = np.percentile(internal_branch_lengths, 75),
            standard_deviation = stat.stdev(internal_branch_lengths),
            variance           = stat.variance(internal_branch_lengths),
            minimum            = np.min(internal_branch_lengths),
            maximum            = np.max(internal_branch_lengths)
        )

        return internal_branch_lengths, stats