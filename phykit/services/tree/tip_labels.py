from typing import Dict

from .base import Tree


class TipLabels(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        tree = self.read_tree_file()
        try:
            print("\n".join([tip.name for tip in tree.get_terminals()]))
        except BrokenPipeError:
            pass

    def process_args(self, args) -> Dict[str, str]:
        return dict(tree_file_path=args.tree)
