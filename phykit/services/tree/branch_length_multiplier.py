import getopt
import logging
import os.path
import sys

from Bio import Phylo

from .base import Tree

class BranchLengthMultiplier(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        tree_mod = self.multiply_branch_lengths_by_factor(tree, self.factor)
        self.write_tree_file(tree_mod, self.output_file_path)

    def process_args(self, args):
        if args.output is None:
            output_file_path = f"{args.tree}.factor_{args.factor}.tre"
        else:
            output_file_path = f"{args.output}"

        return dict(
            tree_file_path=args.tree,
            factor=args.factor,
            output_file_path=output_file_path
        )

    def multiply_branch_lengths_by_factor(self, tree, factor):
        for internode in tree.get_nonterminals():
            if internode.branch_length != None:
                internode.branch_length = internode.branch_length*factor
        for term in tree.get_terminals():
            if term.branch_length != None:
                term.branch_length = term.branch_length*factor
        return tree
