import sys
import getopt
import os.path
import statistics as stat

from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin
import itertools
import numpy as np

from .base import Tree


class RobinsonFouldsDistance(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree_zero = self.read_tree_file()
        tree_one = self.read_tree1_file()
        plain_rf, normalized_rf = self.calculate_robinson_foulds_distance(tree_zero, tree_one)
        
        if (plain_rf, normalized_rf):
            print(f"{plain_rf}\t{normalized_rf}")

    def process_args(self, args):
        return dict(tree_file_path=args.tree_zero, tree1_file_path=args.tree_one)

    def calculate_robinson_foulds_distance(self, tree_zero, tree_one):
        plain_rf = 0
        trees = [[tree_zero, tree_one], [tree_one, tree_zero]]
        for tree in trees:
            # loop through tree_zero and find similar clade in tree_one
            for clade_zero in tree[0].get_nonterminals():
                # initialize and populate a list of tip names in tree_zero
                tip_names_zero = []
                for leaf in clade_zero.get_terminals():
                    tip_names_zero.append(leaf.name)
                # get common ancestor of tree_zero tip names in tree_one
                clade_one = tree[1].common_ancestor(tip_names_zero)
                # initialize and populate a list of tip names in tree_one
                tip_names_one = []
                for leaf in clade_one.get_terminals():
                    tip_names_one.append(leaf.name)

                # compare the list of tip names
                if set(tip_names_zero) == set(tip_names_one):
                    continue
                elif set(tip_names_zero) != set(tip_names_one):
                    plain_rf +=1
        
        # count the number of tips in a phylogeny
        tip_count = 0
        for tip in trees[0][0].get_terminals():
            tip_count += 1

        # calculate normalized rf distance
        normalized_rf = (plain_rf/(2*(tip_count-3)))

        return plain_rf, normalized_rf

