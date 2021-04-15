import getopt
import logging
import math
import os.path
import statistics as stat
import sys

from Bio import Phylo
import numpy as np

from .base import Tree
from ...helpers.stats_summary import calculate_summary_statistics_from_arr, print_summary_statistics

class BipartitionSupportStats(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        bs_vals = self.get_bipartition_support_vals(tree) 

        if self.verbose:
            try:
                for bs_val in bs_vals:
                    print(bs_val)
            except BrokenPipeError:
                pass
        else:
            stats = calculate_summary_statistics_from_arr(bs_vals) 
            print_summary_statistics(stats)

    def process_args(self, args):
        return dict(tree_file_path=args.tree, verbose=args.verbose)

    def get_bipartition_support_vals(self, tree):
        # initialize list to hold bootstrap values
        bs_vals = []

        # populate bs_vals with bootstrap values
        for terminal in tree.get_nonterminals():
            # only include if a bootstrap value is present
            if terminal.confidence != None:
                bs_vals.append(terminal.confidence)
        return bs_vals