from typing import Dict
import pickle

from Bio.Phylo import Newick

from .base import Tree
from ...helpers.json_output import print_json


class BranchLengthMultiplier(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            factor=parsed["factor"],
            output_file_path=parsed["output_file_path"],
        )
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        tree = self.read_tree_file()
        # Make a deep copy to avoid modifying the cached tree
        tree_copy = pickle.loads(pickle.dumps(tree, protocol=pickle.HIGHEST_PROTOCOL))
        scaled_count = self.multiply_branch_lengths_by_factor(tree_copy, self.factor)
        self.write_tree_file(tree_copy, self.output_file_path)

        if self.json_output:
            print_json(
                dict(
                    input_tree=self.tree_file_path,
                    factor=self.factor,
                    scaled_branches=scaled_count,
                    output_file=self.output_file_path,
                )
            )

    def process_args(self, args) -> Dict[str, str]:
        output_file_path = \
            args.output or f"{args.tree}.factor_{args.factor}.tre"
        return dict(
            tree_file_path=args.tree,
            factor=args.factor,
            output_file_path=output_file_path,
            json_output=getattr(args, "json", False),
        )

    def multiply_branch_lengths_by_factor(
        self,
        tree: Newick.Tree,
        factor: float,
    ) -> int:
        scaled_count = 0
        for node in tree.get_nonterminals() + tree.get_terminals():
            if node.branch_length is not None:
                node.branch_length *= factor
                scaled_count += 1
        return scaled_count
