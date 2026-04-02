import pickle
from Bio import Phylo

from .base import Tree

from ...helpers.files import read_single_column_file_to_list
from ...helpers.json_output import print_json


class RootTree(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            outgroup_taxa_file_path=parsed["outgroup_taxa_file_path"],
            output_file_path=parsed["output_file_path"],
        )
        self.json_output = parsed["json_output"]

    def run(self):
        tree = self.read_tree_file()
        # Make a deep copy to avoid modifying the cached tree
        tree_copy = pickle.loads(pickle.dumps(tree, protocol=pickle.HIGHEST_PROTOCOL))

        outgroup = \
            read_single_column_file_to_list(self.outgroup_taxa_file_path)

        Phylo.BaseTree.Tree.root_with_outgroup(tree_copy, outgroup)

        self.write_tree_file(tree_copy, self.output_file_path)

        if self.json_output:
            print_json(
                dict(
                    input_tree=self.tree_file_path,
                    outgroup_taxa_file=self.outgroup_taxa_file_path,
                    outgroup_taxa=outgroup,
                    output_file=self.output_file_path,
                )
            )

    def process_args(self, args):
        tree_file_path = args.tree

        output_file_path = \
            args.output if args.output else f"{tree_file_path}.rooted"

        return dict(
            tree_file_path=tree_file_path,
            outgroup_taxa_file_path=args.root,
            output_file_path=output_file_path,
            json_output=getattr(args, "json", False),
        )
