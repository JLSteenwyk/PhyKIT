import getopt
import logging
import math
import os.path
import statistics as stat
import sys

from Bio import Phylo
import numpy as np


from phykit.services.tree.base import Tree

class InternalBranchStats(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        mean, median, twenty_fifth, seventy_fifth, minimum, maximum, standard_deviation, variance = self.calculate_internal_branch_stats(tree)
        if (mean, median, twenty_fifth, seventy_fifth, standard_deviation, variance):
            print(f"mean: {mean}")
            print(f"median: {median}")
            print(f"25th percentile: {twenty_fifth}")
            print(f"75th percentile: {seventy_fifth}")
            print(f"minimum: {minimum}")
            print(f"maximum: {maximum}")
            print(f"standard deviation: {standard_deviation}")
            print(f"variance: {variance}")

    def process_args(self, args):
        return dict(tree_file_path=args.tree)

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
        
        mean               = stat.mean(internal_branch_lengths)
        median             = stat.median(internal_branch_lengths)
        twenty_fifth       = np.percentile(internal_branch_lengths, 25)
        seventy_fifth      = np.percentile(internal_branch_lengths, 75)
        standard_deviation = stat.stdev(internal_branch_lengths)
        variance           = stat.variance(internal_branch_lengths)
        minimum            = np.min(internal_branch_lengths)
        maximum            = np.max(internal_branch_lengths)

        return mean, median, twenty_fifth, seventy_fifth, minimum, maximum, standard_deviation, variance