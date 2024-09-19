from typing import Dict

from Bio.Phylo import Newick

from .base import Tree


class BranchLengthMultiplier(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        tree = self.read_tree_file()
        self.multiply_branch_lengths_by_factor(tree, self.factor)
        self.write_tree_file(tree, self.output_file_path)

    def process_args(self, args) -> Dict[str, str]:
        output_file_path = \
            args.output or f"{args.tree}.factor_{args.factor}.tre"
        return dict(
            tree_file_path=args.tree,
            factor=args.factor,
            output_file_path=output_file_path,
        )

    def multiply_branch_lengths_by_factor(
        self,
        tree: Newick.Tree,
        factor: float,
    ) -> Newick.Tree:
        for node in tree.get_nonterminals() + tree.get_terminals():
            if node.branch_length is not None:
                node.branch_length *= factor
        return tree
