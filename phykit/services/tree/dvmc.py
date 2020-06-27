import getopt
import logging
import math
import os.path
import statistics as stat
import sys

from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin
import itertools
import numpy as np

from .base import Tree

from ...helpers.files import read_single_column_file_to_list

class DVMC(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        outgroup = read_single_column_file_to_list(self.outgroup_taxa_file_path)
        dvmc = self.calculate_dvmc(tree, outgroup)

        print(f"{dvmc}")

    def process_args(self, args):
        return dict(tree_file_path=args.tree, outgroup_taxa_file_path=args.root)

    def get_names_of_outgroup_taxa_that_are_present(
        self,
        outgroup: list,
        tree: Tree
    ) -> list:
        # initialize list for outgroup taxa that are present in tree
        out_pres = []
        # loop through terminal branch
        for term in tree.get_terminals():
            # if outgroup taxa is present, append it to out_pres
            if term.name in outgroup:
                out_pres.append(term.name)

        return out_pres

    def calculate_dvmc(self, tree: Tree, outgroup: list):
        # get names of outgroup taxa in tree
        out_pres = self.get_names_of_outgroup_taxa_that_are_present(outgroup, tree)

        # prune outgroup taxa from tree
        tree = self.prune_tree_using_taxa_list(tree, outgroup)

        # determine number of taxa in the tree
        num_spp = tree.count_terminals()

        # initialize list of root to tip distances, the sum
        # of distances in the tree, the average of distances
        # in the tree
        dist = []
        sum_dist = float(0.0)
        avg_dist = float(0.0)

        # loop through terminal branches and store distances from
        # the root to the tip in a list
        for term in tree.get_terminals():
            # append root to tip distance to dist
            dist.append((tree.distance(term)))
            # keep running sum of tree distances
            sum_dist += tree.distance(term)

        # calculate average tree distance
        avg_dist = sum_dist/num_spp

        # determine the sum of i=1 to N for (x_i-x_bar)^2
        sumi2N = float(0.0)
        for x_i in dist:
            sumi2N += ((x_i-avg_dist)**2)

        # multiply sumi2N by 1/(N-1) where N is the number of spp
        # and then take the square root
        dvmc = float(0.0)
        dvmc = math.sqrt((1/(num_spp-1))*sumi2N) 

        # return result
        return dvmc