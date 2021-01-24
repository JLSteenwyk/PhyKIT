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
        
        try:
            print(f"{plain_rf}\t{round(normalized_rf, 4)}")
        except BrokenPipeError:
            pass

    def process_args(self, args):
        return dict(tree_file_path=args.tree_zero, tree1_file_path=args.tree_one)

    def calculate_robinson_foulds_distance(self, tree_zero, tree_one):
        plain_rf = 0
        plain_rf = self.compare_trees(plain_rf, tree_zero, tree_one)
        plain_rf = self.compare_trees(plain_rf, tree_one, tree_zero)
        # count the number of tips in a phylogeny
        tip_count = tree_zero.count_terminals()
        
        # calculate normalized rf distance
        normalized_rf = (plain_rf/(2*(tip_count-3)))

        return plain_rf, normalized_rf

    def compare_trees(
        self,
        plain_rf: int,
        tree_zero: Tree,
        tree_one: Tree
    ) -> int:
        # loop through tree_zero and find similar clade in tree_one
        for clade_zero in tree_zero.get_nonterminals():
            # initialize and populate a list of tip names in tree_zero
            tip_names_zero = self.get_tip_names_from_tree(clade_zero)
            # get common ancestor of tree_zero tip names in tree_one
            clade_one = tree_one.common_ancestor(tip_names_zero)
            # initialize and populate a list of tip names in tree_one
            tip_names_one = self.get_tip_names_from_tree(clade_one)

            # compare the list of tip names
            plain_rf = self.determine_if_clade_differs(
                plain_rf,
                tip_names_zero,
                tip_names_one
            )

        return plain_rf

    def determine_if_clade_differs(
        self,
        plain_rf: int,
        tip_names_zero: list,
        tip_names_one: list
    ) -> int:
        """
        if clade differs, add 1 to plain_rf value
        """
        if set(tip_names_zero) != set(tip_names_one):
            plain_rf +=1
        
        return plain_rf



