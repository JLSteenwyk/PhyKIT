import sys
import getopt
import os.path
import statistics as stat
import logging

from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin
import itertools
import numpy as np

from phykit.services.tree.base import Tree

# TODO: how do I read in multiple args?
class DVMC(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        print(tree)
        #outgroup = [line.rstrip('\n') for line in open(outgroup_taxa_file_path)]

        mean, median, twenty_fifth, seventy_fifth, standard_deviation, variance = self.calculate_lb_score(tree)
        if (mean, median, twenty_fifth, seventy_fifth, standard_deviation, variance):
            print(f"mean: {mean}")
            print(f"median: {median}")
            print(f"25th percentile: {twenty_fifth}")
            print(f"75th percentile: {seventy_fifth}")
            print(f"standard deviation: {standard_deviation}")
            print(f"variance: {variance}")

    def process_args(self, args):
        root_taxa_file = args.root
        if not os.path.isfile(root_taxa_file):
            logger.warning("Input file of taxa to root the tree on does not exist.")
            sys.exit()
        return dict(tree_file_path=args.tree, outgroup_taxa_file_path=root_taxa_file)

    def calculate_dvmc(self, tree):
        print("entering calculating dvmc function")
        # # get tree tips
        # tips = []
        # for tip in tree.get_terminals():
        #     tips.append(tip.name)
        
        # # determine pairwise combinations of tips
        # combos = list(itertools.combinations(tips, 2))

        # # determine average distance between tips
        # # avg_dist is PDa
        # avg_dist = float()
        # for combo in combos:
        #     avg_dist+=tree.distance(combo[0], combo[1])
        # avg_dist = avg_dist/len(combos)

        # # calculate average distance of taxon i to all other taxa
        # # or PDi and save each result to LBi
        # avg_PDis = []
        # for tip in tips:
        #     tips_minus_i = list(set(tips) - set(tip))
        #     PDi = []
        #     for tip_minus in tips_minus_i:    
        #         PDi.append(tree.distance(tip, tip_minus))
        #     PDi = sum(PDi) / len(PDi)
        #     avg_PDis.append(PDi)

        # # use PDis and avgDist to calculate LB values for each taxon
        # LBis = []
        # for PDi in avg_PDis:
        #     try:
        #         LBis.append((((PDi/avg_dist)-1)*100))
        #     except ZeroDivisionError:
        #         print("Invalid tree. Tree should contain branch lengths")
        #         return None

        # ## TODO: fix code to write out output of lb scores per taxa
        # ## output should have the same name as the input file
        # # with open(tree + ".LBi-scores.txt", 'w') as f:
        # #     for tip, LBi in zip(tips, LBis):
        # #         f.write(str(tip) + "\t" + str(LBi) + "\n")

        # mean               = stat.mean(LBis)
        # median             = stat.median(LBis)
        # twenty_fifth       = np.percentile(LBis, 25)
        # seventy_fifth      = np.percentile(LBis, 75)
        # standard_deviation = stat.stdev(LBis)
        # variance           = stat.variance(LBis)


        # return mean, median, twenty_fifth, seventy_fifth, standard_deviation, variance