import getopt
import logging
import math
import os.path
import statistics as stat
import sys

from Bio import Phylo
import numpy as np


from phykit.services.tree.base import Tree

class BipartitionSupportStats(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        mean, median, twenty_fifth, seventy_fifth, minimum, maximum, standard_deviation, variance = self.calculate_bipartition_support_stats(tree)
        if (mean, median, twenty_fifth, seventy_fifth, minimum, maximum, standard_deviation, variance):
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

    def calculate_bipartition_support_stats(self, tree):
        # initialize list to hold bootstrap values
        bs_vals = []

        # populate bs_vals with bootstrap values
        for terminal in tree.get_nonterminals():
            # only include if a bootstrap value is present
            if terminal.confidence != None:
                bs_vals.append(terminal.confidence)
        
        mean               = stat.mean(bs_vals)
        median             = stat.median(bs_vals)
        twenty_fifth       = np.percentile(bs_vals, 25)
        seventy_fifth      = np.percentile(bs_vals, 75)
        standard_deviation = stat.stdev(bs_vals)
        variance           = stat.variance(bs_vals)
        minimum            = np.min(bs_vals)
        maximum            = np.max(bs_vals)

        return mean, median, twenty_fifth, seventy_fifth, minimum, maximum, standard_deviation, variance