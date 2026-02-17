from typing import Dict
import copy

from .base import Tree
from ...helpers.json_output import print_json


class CollapseBranches(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            support=parsed["support"],
            output_file_path=parsed["output_file_path"],
        )
        self.json_output = parsed["json_output"]

    def run(self):
        tree = self.read_tree_file()
        # Make a deep copy to avoid modifying the cached tree
        tree_copy = copy.deepcopy(tree)
        initial_nonterminals = len(tree_copy.get_nonterminals())
        tree_copy.collapse_all(
            lambda c: c.confidence and c.confidence < self.support
        )
        final_nonterminals = len(tree_copy.get_nonterminals())
        self.write_tree_file(tree_copy, self.output_file_path)

        if self.json_output:
            print_json(
                dict(
                    input_tree=self.tree_file_path,
                    support_threshold=self.support,
                    collapsed_branches=max(0, initial_nonterminals - final_nonterminals),
                    output_file=self.output_file_path,
                )
            )

    def process_args(self, args) -> Dict[str, str]:
        output_file_path = \
            args.output or f"{args.tree}.collapsed_{args.support}.tre"
        return dict(
            tree_file_path=args.tree,
            support=args.support,
            output_file_path=output_file_path,
            json_output=getattr(args, "json", False),
        )
