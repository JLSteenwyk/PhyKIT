import logging

from Bio import Phylo

from .base import Tree

from ...helpers.files import read_single_column_file_to_list

class RootTree(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        
        outgroup = read_single_column_file_to_list(self.outgroup_taxa_file_path)

        Phylo.BaseTree.Tree.root_with_outgroup(tree, outgroup)

        self.write_tree_file(tree, self.output_file_path)
    
    def process_args(self, args):
        tree_file_path=args.tree

        if args.output is None:
            output_file_path = f"{tree_file_path}.rooted"
        else:
            output_file_path = f"{args.output}"

        return dict(
            tree_file_path=tree_file_path,
            outgroup_taxa_file_path=args.root,
            output_file_path=output_file_path
        )