import sys
import getopt
import os.path
import statistics as stat

from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin
import itertools
import numpy as np

from .base import Tree


class PatristicDistances(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        mean, median, twenty_fifth, seventy_fifth, minimum, maximum, standard_deviation, variance, patristic_distances, combos = self.calculate_patristic_distances(tree)
        if not self.verbose:
            if (mean, median, twenty_fifth, seventy_fifth, minimum, maximum, standard_deviation, variance):
                print(f"mean: {mean}")
                print(f"median: {median}")
                print(f"25th percentile: {twenty_fifth}")
                print(f"75th percentile: {seventy_fifth}")
                print(f"minimum: {minimum}")
                print(f"maximum: {maximum}")
                print(f"standard deviation: {standard_deviation}")
                print(f"variance: {variance}")
        elif self.verbose:
            for combo, patristic_distance in zip(combos, patristic_distances):
                print(f"{combo[0]}-{combo[1]}\t{patristic_distance}")

    def process_args(self, args):
        return dict(tree_file_path=args.tree, verbose=args.verbose)

    def calculate_patristic_distances(self, tree):
        # get tree tips
        tips = []
        for tip in tree.get_terminals():
            tips.append(tip.name)
        
        # determine pairwise combinations of tips
        combos = list(itertools.combinations(tips, 2))

        # determine average distance between tips
        patristic_distances = []
        for combo in combos:
            patristic_distances.append(tree.distance(combo[0], combo[1]))

        mean               = stat.mean(patristic_distances)
        median             = stat.median(patristic_distances)
        twenty_fifth       = np.percentile(patristic_distances, 25)
        seventy_fifth      = np.percentile(patristic_distances, 75)
        minimum            = np.min(patristic_distances)
        maximum            = np.max(patristic_distances)
        standard_deviation = stat.stdev(patristic_distances)
        variance           = stat.variance(patristic_distances)


        return mean, median, twenty_fifth, seventy_fifth, minimum, maximum, standard_deviation, variance, patristic_distances, combos