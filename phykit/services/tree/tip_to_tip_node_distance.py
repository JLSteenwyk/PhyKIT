import sys
from typing import Dict

from Bio.Phylo.BaseTree import TreeMixin
from Bio.Phylo import Newick

from .base import Tree
from ...helpers.json_output import print_json


class TipToTipNodeDistance(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            tip_1=parsed["tip_1"],
            tip_2=parsed["tip_2"],
        )
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        tree_zero = self.read_tree_file()

        self.check_leaves(tree_zero, self.tip_1, self.tip_2)

        node_distance = len(TreeMixin.trace(tree_zero, self.tip_1, self.tip_2))
        if self.json_output:
            print_json(
                dict(
                    taxon_a=self.tip_1,
                    taxon_b=self.tip_2,
                    tip_to_tip_node_distance=node_distance,
                )
            )
            return
        print(node_distance)

    def check_leaves(
        self,
        tree_zero: Newick.Tree,
        tip_1: str,
        tip_2: str,
    ) -> None:
        leaf1 = TreeMixin.find_any(tree_zero, tip_1)
        if not bool(leaf1):
            print(tip_1, "not on tree\nExiting...")
            sys.exit(2)
        leaf2 = TreeMixin.find_any(tree_zero, tip_2)
        if not bool(leaf2):
            print(tip_2, "not on tree\nExiting...")
            sys.exit(2)

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            tree_file_path=args.tree_zero,
            tip_1=args.tip_1,
            tip_2=args.tip_2,
            json_output=getattr(args, "json", False),
        )
