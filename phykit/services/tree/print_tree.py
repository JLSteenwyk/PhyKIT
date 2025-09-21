from typing import Dict
import copy

from Bio import Phylo

from .base import Tree


class PrintTree(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()

        if self.remove:
            # Make a deep copy to avoid modifying the cached tree
            tree = copy.deepcopy(tree)
            for node in tree.get_terminals() + tree.get_nonterminals():
                node.branch_length = None

        try:
            Phylo.draw_ascii(tree)
        except BrokenPipeError:
            pass

    def process_args(self, args) -> Dict[str, str]:
        return dict(tree_file_path=args.tree, remove=args.remove)
