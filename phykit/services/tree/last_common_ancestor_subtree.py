import sys

from Bio import Phylo

from .base import Tree

from ...helpers.files import read_single_column_file_to_list


class LastCommonAncestorSubtree(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        try:
            taxa = read_single_column_file_to_list(self.list_of_taxa)
        except FileNotFoundError:
            print("Taxa list file is not found. Please check pathing.")
            sys.exit()
        tree = tree.common_ancestor(taxa)

        self.write_tree_file(tree, self.output_file_path)

    def process_args(self, args):
        tree_file_path = args.tree

        if args.output is None:
            output_file_path = f"{tree_file_path}.subtree.tre"
        else:
            output_file_path = f"{args.output}"

        return dict(
            tree_file_path=tree_file_path,
            output_file_path=output_file_path,
            list_of_taxa=args.list_of_taxa,
        )
