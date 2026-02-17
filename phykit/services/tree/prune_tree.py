from typing import Dict
import copy

from .base import Tree

from ...helpers.files import read_single_column_file_to_list
from ...helpers.json_output import print_json


class PruneTree(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            list_of_taxa=parsed["list_of_taxa"],
            output_file_path=parsed["output_file_path"],
            keep=parsed["keep"],
        )
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        tree = self.read_tree_file()
        # Make a deep copy to avoid modifying the cached tree
        tree_copy = copy.deepcopy(tree)

        taxa = read_single_column_file_to_list(self.list_of_taxa)

        if self.keep:
            tips_in_tree = [term.name for term in tree_copy.get_terminals()]
            taxa = [x for x in tips_in_tree if x not in taxa]

        tree_copy = self.prune_tree_using_taxa_list(tree_copy, taxa)

        self.write_tree_file(tree_copy, self.output_file_path)

        if self.json_output:
            print_json(
                dict(
                    input_tree=self.tree_file_path,
                    input_taxa_file=self.list_of_taxa,
                    keep_input_taxa=self.keep,
                    taxa_pruned=sorted(taxa),
                    pruned_count=len(taxa),
                    remaining_tips=tree_copy.count_terminals(),
                    output_file=self.output_file_path,
                )
            )

    def process_args(self, args) -> Dict[str, str]:
        tree_file_path = args.tree
        output_file_path = \
            f"{args.output}" if args.output else f"{tree_file_path}.pruned"

        keep = True if args.keep is None else args.keep

        return dict(
            tree_file_path=tree_file_path,
            list_of_taxa=args.list_of_taxa,
            output_file_path=output_file_path,
            keep=keep,
            json_output=getattr(args, "json", False),
        )
