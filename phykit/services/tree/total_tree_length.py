from typing import Dict, Union

from Bio.Phylo import Newick
from .base import Tree
from ...helpers.json_output import print_json


class TotalTreeLength(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        tree = self.read_tree_file()
        total_tree_length = round(self.calculate_total_tree_length(tree), 4)
        if self.json_output:
            print_json(dict(total_tree_length=total_tree_length))
            return
        print(total_tree_length)

    def process_args(self, args) -> Dict[str, str]:
        return dict(tree_file_path=args.tree, json_output=getattr(args, "json", False))

    def calculate_total_tree_length(
        self,
        tree: Newick.Tree
    ) -> Union[int, float]:
        total_len = tree.total_branch_length()

        if isinstance(total_len, (int, float)):
            return total_len
