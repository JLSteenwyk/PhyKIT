import copy
from Bio import Phylo

from .base import Tree

from ...helpers.files import read_single_column_file_to_list


class RootTree(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        # Make a deep copy to avoid modifying the cached tree
        tree_copy = copy.deepcopy(tree)

        outgroup = \
            read_single_column_file_to_list(self.outgroup_taxa_file_path)

        Phylo.BaseTree.Tree.root_with_outgroup(tree_copy, outgroup)

        self.write_tree_file(tree_copy, self.output_file_path)

    def process_args(self, args):
        tree_file_path = args.tree

        output_file_path = \
            args.output if args.output else f"{tree_file_path}.rooted"

        return dict(
            tree_file_path=tree_file_path,
            outgroup_taxa_file_path=args.root,
            output_file_path=output_file_path,
        )
