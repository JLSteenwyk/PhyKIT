import sys
import getopt
import os.path
import statistics as stat

from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin
import itertools
import numpy as np

from .base import Tree

from ...helpers.stats_summary import calculate_summary_statistics_from_arr, print_summary_statistics


class LBScore(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        tips, LBis = self.calculate_lb_score(tree)
        if self.verbose:
            for tip, LBi in zip(tips, LBis):
                print(f"{tip}\t{LBi}")
        else:
            stats = calculate_summary_statistics_from_arr(LBis)
            print_summary_statistics(stats)

    def process_args(self, args):
        return dict(tree_file_path=args.tree, verbose=args.verbose)

    def calculate_lb_score(self, tree):
        # get tree tips
        tips = []
        for tip in tree.get_terminals():
            tips.append(tip.name)
        
        # determine pairwise combinations of tips
        combos = list(itertools.combinations(tips, 2))

        # determine average distance between tips
        # avg_dist is PDa
        avg_dist = float()
        for combo in combos:
            avg_dist+=tree.distance(combo[0], combo[1])
        avg_dist = avg_dist/len(combos)

        # calculate average distance of taxon i to all other taxa
        # or PDi and save each result to LBi
        avg_PDis = []
        for tip in tips:
            tips_minus_i = list(set(tips) - set(tip))
            PDi = []
            for tip_minus in tips_minus_i:    
                PDi.append(tree.distance(tip, tip_minus))
            PDi = sum(PDi) / len(PDi)
            avg_PDis.append(PDi)

        # use PDis and avgDist to calculate LB values for each taxon
        LBis = []
        for PDi in avg_PDis:
            try:
                LBis.append((((PDi/avg_dist)-1)*100))
            except ZeroDivisionError:
                print("Invalid tree. Tree should contain branch lengths")
                return None
        
        return tips, LBis
