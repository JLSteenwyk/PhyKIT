import logging

from Bio import Phylo

from .base import Tree

class PrintTree(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        
        # remove branch lengths if specified by user
        # otherwise, print ascii tree
        if self.remove:
            for internode in tree.get_nonterminals():
                internode.branch_length = None
            for leaf in tree.get_terminals():
                leaf.branch_length = None

        try:
            Phylo.draw_ascii(tree)
        except BrokenPipeError:
            pass
    
    def process_args(self, args):
        return dict(tree_file_path=args.tree, remove=args.remove)