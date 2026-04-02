from typing import Dict
import pickle

from Bio.Phylo import Newick

from .base import Tree
from ...helpers.json_output import print_json


class InternodeLabeler(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            output_file_path=parsed["output_file_path"],
        )
        self.json_output = parsed["json_output"]

    def run(self):
        tree = self.read_tree_file()
        # Make a deep copy to avoid modifying the cached tree
        tree_copy = pickle.loads(pickle.dumps(tree, protocol=pickle.HIGHEST_PROTOCOL))
        labeled_count = self.add_labels_to_tree(tree_copy)
        self.write_tree_file(tree_copy, self.output_file_path)

        if self.json_output:
            print_json(
                dict(
                    input_tree=self.tree_file_path,
                    labeled_internodes=labeled_count,
                    output_file=self.output_file_path,
                )
            )

    def process_args(self, args) -> Dict[str, str]:
        output_file_path = args.output or f"{args.tree}.internode_labels.tre"

        return dict(
            tree_file_path=args.tree,
            output_file_path=output_file_path,
            json_output=getattr(args, "json", False),
        )

    def add_labels_to_tree(
        self,
        tree: Newick.Tree
    ) -> int:
        labeled_count = 0
        for label, node in enumerate(tree.get_nonterminals(), start=1):
            node.confidence = label
            labeled_count += 1
        return labeled_count
