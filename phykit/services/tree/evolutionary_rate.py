from typing import Dict

from .base import Tree


class EvolutionaryRate(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        tree = self.read_tree_file()
        total_tree_length = tree.total_branch_length()
        num_terminals = tree.count_terminals()
        print(round(total_tree_length / num_terminals, 4))

    def process_args(self, args) -> Dict[str, str]:
        return dict(tree_file_path=args.tree)
