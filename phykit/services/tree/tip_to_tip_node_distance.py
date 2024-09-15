import sys
from typing import Dict

from Bio.Phylo.BaseTree import TreeMixin
from Bio.Phylo import Newick

from .base import Tree


class TipToTipNodeDistance(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        tree_zero = self.read_tree_file()

        self.check_leaves(tree_zero, self.tip_1, self.tip_2)

        print(len(TreeMixin.trace(tree_zero, self.tip_1, self.tip_2)))

    def check_leaves(
        self,
        tree_zero: Newick.Tree,
        tip_1: str,
        tip_2: str,
    ) -> None:
        leaf1 = TreeMixin.find_any(tree_zero, tip_1)
        if not bool(leaf1):
            print(tip_1, "not on tree\nExiting...")
            sys.exit()
        leaf2 = TreeMixin.find_any(tree_zero, tip_2)
        if not bool(leaf2):
            print(tip_2, "not on tree\nExiting...")
            sys.exit()

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            tree_file_path=args.tree_zero,
            tip_1=args.tip_1,
            tip_2=args.tip_2,
        )
