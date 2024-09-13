from typing import Dict

from .base import Tree


class CollapseBranches(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        tree.collapse_all(
            lambda c: c.confidence and c.confidence < self.support
        )
        self.write_tree_file(tree, self.output_file_path)

    def process_args(self, args) -> Dict[str, str]:
        output_file_path = \
            args.output or f"{args.tree}.collapsed_{args.support}.tre"
        return dict(
            tree_file_path=args.tree,
            support=args.support,
            output_file_path=output_file_path,
        )
