import logging

from Bio import Phylo

from .base import Tree

class PruneTree(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        
        taxa = [line.rstrip('\n') for line in open(self.list_of_taxa)]

        for taxon in taxa:
            tree.prune(taxon)

        self.write_tree_file(tree, self.output_file_path)
    
    def process_args(self, args):
        return dict(
            tree_file_path=args.tree,
            list_of_taxa=args.list_of_taxa,
            output_file_path=f"{args.tree}.pruned",
        )