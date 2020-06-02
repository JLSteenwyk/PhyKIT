from Bio import Phylo

from ..base import BaseService


class Tree(BaseService):
    def __init__(
        self,
        *args,
        tree_file_path,
        outgroup_taxa_file_path=None,
        output_file_path=None,
        verbose=None
    ):
        self.tree_file_path = tree_file_path
        self.output_file_path = output_file_path
        self.outgroup_taxa_file_path = outgroup_taxa_file_path
        self.tree_format = "newick"
        self.verbose = verbose

    def read_tree_file(self):
        return Phylo.read(self.tree_file_path, self.tree_format)

    def write_tree_file(self, tree, output_file_path):
        return Phylo.write(tree, output_file_path, self.tree_format)
