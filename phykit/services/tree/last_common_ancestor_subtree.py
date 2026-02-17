import sys
import copy
from typing import Dict

from .base import Tree

from ...helpers.files import read_single_column_file_to_list
from ...helpers.json_output import print_json


class LastCommonAncestorSubtree(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            output_file_path=parsed["output_file_path"],
            list_of_taxa=parsed["list_of_taxa"],
        )
        self.json_output = parsed["json_output"]

    def run(self):
        tree = self.read_tree_file()
        # Make a deep copy to avoid issues with cached tree modifications
        tree_copy = copy.deepcopy(tree)
        try:
            taxa = read_single_column_file_to_list(self.list_of_taxa)
        except FileNotFoundError:
            print("Taxa list file is not found. Please check pathing.")
            sys.exit(2)
        subtree = tree_copy.common_ancestor(taxa)

        self.write_tree_file(subtree, self.output_file_path)

        if self.json_output:
            print_json(
                dict(
                    input_tree=self.tree_file_path,
                    input_taxa_file=self.list_of_taxa,
                    taxa_count=len(taxa),
                    subtree_tips=subtree.count_terminals(),
                    output_file=self.output_file_path,
                )
            )

    def process_args(self, args) -> Dict[str, str]:
        tree_file_path = args.tree
        output_file_path = args.output or f"{tree_file_path}.subtree.tre"

        return dict(
            tree_file_path=tree_file_path,
            output_file_path=output_file_path,
            list_of_taxa=args.list_of_taxa,
            json_output=getattr(args, "json", False),
        )
