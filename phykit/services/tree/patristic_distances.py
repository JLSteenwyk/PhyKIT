import sys
import getopt
import os.path
import statistics as stat
from typing import Tuple

from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin
import itertools
import numpy as np

from .base import Tree

from ...helpers.stats_summary import calculate_summary_statistics_from_arr, print_summary_statistics

class PatristicDistances(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        patristic_distances, combos, stats = self.calculate_patristic_distances(tree)
        if self.verbose:
            try:
                for combo, patristic_distance in zip(combos, patristic_distances):
                    print(f"{combo[0]}\t{combo[1]}\t{round(patristic_distance, 4)}")
            except BrokenPipeError:
                pass
        else:
            print_summary_statistics(stats)

    def process_args(self, args):
        return dict(tree_file_path=args.tree, verbose=args.verbose)

    def calculate_distance_between_pairs(self, tips: list, tree) -> Tuple[list, list]:
        # determine pairwise combinations of tips
        combos = list(itertools.combinations(tips, 2))

        # determine average distance between tips
        patristic_distances = []
        for combo in combos:
            patristic_distances.append(tree.distance(combo[0], combo[1]))

        return combos, patristic_distances

    def calculate_patristic_distances(self, tree):
        # get tree tips
        tips = self.get_tip_names_from_tree(tree)
        
        # get distances between pairs of taxa
        combos, patristic_distances = self.calculate_distance_between_pairs(tips, tree)
        
        # calculate summary stats
        stats = calculate_summary_statistics_from_arr(patristic_distances)
        
        return patristic_distances, combos, stats