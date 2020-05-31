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

class DVMC(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        outgroup = [line.rstrip('\n') for line in open(self.outgroup_taxa_file_path)]
        dvmc = self.calculate_dvmc(tree, outgroup)

        if dvmc:
            print(f"DVMC: {dvmc}")

    def process_args(self, args):
        return dict(tree_file_path=args.tree, outgroup_taxa_file_path=args.root)

    def calculate_dvmc(self, tree, outgroup):
        # initialize list for outgroup taxa present in the tree
        out_pres   = []

        # initialize value to hold the number of species
        num_spp = float(0.0)

        # loop through terminal branch
        for term in tree.get_terminals():
            # if outgroup taxa is present, append it to out_pres
            if term.name in outgroup:
                out_pres.append(term.name)

        # root tree on outgroup
        tree.root_with_outgroup(out_pres)

        # prune outgroiup taxa from the tree
        for taxon in out_pres:
            tree.prune(taxon)

        # determine number of taxa in the tree
        for term in tree.get_terminals():
            # add 1 for each tip that is iterated through
            num_spp+=1

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

        # multiple sumi2N by 1/(N-1) where N is the number of spp
        # and take the square root
        dvmc = float(0.0)
        dvmc = math.sqrt((1/(num_spp-1))*sumi2N) 

        # return result
        return dvmc