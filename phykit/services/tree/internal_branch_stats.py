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

class InternalBranchStats(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        internal_branch_lengths, stats = self.calculate_internal_branch_stats(tree)

        if self.verbose:
            try:
                for internal_branch_length in internal_branch_lengths:
                    print(round(internal_branch_length, 4))
            except BrokenPipeError:
                pass
        else:
            print_summary_statistics(stats)

    def process_args(self, args):
        return dict(tree_file_path=args.tree, verbose=args.verbose)

    def get_internal_branch_lengths(self, tree) -> list:
        """
        loop through tree and get all internal branch lengths
        """
        internal_branch_lengths = []
        for internal_branch in tree.get_nonterminals():
            if internal_branch.branch_length != None:
                internal_branch_lengths.append(internal_branch.branch_length)

        return internal_branch_lengths

    def check_tree_has_branch_lengths(self, internal_branch_lengths:list) -> None:
        """
        if tree has no branch lengths, exit
        """
        if len(internal_branch_lengths) == 0:
            print("Calculating internal branch statistics requires a phylogeny with branch lengths.")
            sys.exit()


    def calculate_internal_branch_stats(self, tree):
        # save internal branch lengths to internal_branch_lengths
        internal_branch_lengths = self.get_internal_branch_lengths(tree)
        
        # If the phylogeny had no branch lengths, inform user and quit
        self.check_tree_has_branch_lengths(internal_branch_lengths)
        
        # calculate summary stats
        stats = calculate_summary_statistics_from_arr(internal_branch_lengths)

        return internal_branch_lengths, stats