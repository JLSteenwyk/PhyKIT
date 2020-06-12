import getopt
import logging
import math
import os.path
import statistics as stat
import sys

from Bio import Phylo
import numpy as np

from .base import Tree

class BipartitionSupportStats(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        bs_vals, stats = self.calculate_bipartition_support_stats(tree)
        if self.verbose:
            for bs_val in bs_vals:
                print(bs_val)
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

    def calculate_bipartition_support_stats(self, tree):
        # initialize list to hold bootstrap values
        bs_vals = []

        # populate bs_vals with bootstrap values
        for terminal in tree.get_nonterminals():
            # only include if a bootstrap value is present
            if terminal.confidence != None:
                bs_vals.append(terminal.confidence)
        
        stats = dict(
            mean               = stat.mean(bs_vals),
            median             = stat.median(bs_vals),
            twenty_fifth       = np.percentile(bs_vals, 25),
            seventy_fifth      = np.percentile(bs_vals, 75),
            standard_deviation = stat.stdev(bs_vals),
            variance           = stat.variance(bs_vals),
            minimum            = np.min(bs_vals),
            maximum            = np.max(bs_vals)
        )

        return bs_vals, stats