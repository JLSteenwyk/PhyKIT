from typing import Dict

from .base import Tree
from ...helpers.json_output import print_json


class Treeness(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        tree = self.read_tree_file()
        treeness = self.calculate_treeness(tree)
        treeness = round(treeness, 4)
        if self.json_output:
            print_json(dict(treeness=treeness))
            return
        print(treeness)

    def process_args(self, args) -> Dict[str, str]:
        return dict(tree_file_path=args.tree, json_output=getattr(args, "json", False))
