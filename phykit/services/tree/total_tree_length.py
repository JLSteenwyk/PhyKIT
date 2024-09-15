from typing import Dict, Union

from Bio.Phylo import Newick
from .base import Tree


class TotalTreeLength(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        tree = self.read_tree_file()
        print(round(self.calculate_total_tree_length(tree), 4))

    def process_args(self, args) -> Dict[str, str]:
        return dict(tree_file_path=args.tree)

    def calculate_total_tree_length(
        self,
        tree: Newick.Tree
    ) -> Union[int, float]:
        total_len = tree.total_branch_length()

        if isinstance(total_len, (int, float)):
            return total_len
