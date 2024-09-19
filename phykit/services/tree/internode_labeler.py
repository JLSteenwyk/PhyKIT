from typing import Dict

from Bio.Phylo import Newick

from .base import Tree


class InternodeLabeler(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        self.add_labels_to_tree(tree)  
        self.write_tree_file(tree, self.output_file_path)

    def process_args(self, args) -> Dict[str, str]:
        output_file_path = args.output or f"{args.tree}.internode_labels.tre"

        return dict(
            tree_file_path=args.tree,
            output_file_path=output_file_path,
        )

    def add_labels_to_tree(
        self,
        tree: Newick.Tree
    ) -> None:
        for label, node in enumerate(tree.get_nonterminals(), start=1):
            node.confidence = label
