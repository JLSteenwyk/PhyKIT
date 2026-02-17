from typing import Dict

from .base import Tree
from ...helpers.json_output import print_json


class EvolutionaryRate(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        tree = self.read_tree_file()
        total_tree_length = tree.total_branch_length()
        num_terminals = tree.count_terminals()
        evo_rate = round(total_tree_length / num_terminals, 4)
        if self.json_output:
            print_json(dict(evolutionary_rate=evo_rate))
            return
        print(evo_rate)

    def process_args(self, args) -> Dict[str, str]:
        return dict(tree_file_path=args.tree, json_output=getattr(args, "json", False))
