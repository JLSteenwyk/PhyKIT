from typing import Dict

from .base import Tree
from ...helpers.json_output import print_json


class TipLabels(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        tree = self.read_tree_file()
        tips = [tip.name for tip in tree.get_terminals()]
        if self.json_output:
            rows = [dict(taxon=tip) for tip in tips]
            print_json(dict(rows=rows, tips=tips))
            return
        try:
            print("\n".join(tips))
        except BrokenPipeError:
            pass

    def process_args(self, args) -> Dict[str, str]:
        return dict(tree_file_path=args.tree, json_output=getattr(args, "json", False))
