from Bio import Phylo

from ..base import BaseService

class Tree(BaseService):
    def __init__(
        self,
        *args,
        tree_file_path,
        tree=None,
        alignment_file_path=None,
        tree1_file_path=None,
        outgroup_taxa_file_path=None,
        output_file_path=None,
        factor=None,
        verbose=None
    ):
        self.tree_file_path = tree_file_path
        self.tree1_file_path = tree1_file_path
        self.tree = tree
        self.alignment_file_path = alignment_file_path
        self.output_file_path = output_file_path
        self.outgroup_taxa_file_path = outgroup_taxa_file_path
        self.tree_format = "newick"
        self.verbose = verbose
        self.factor = factor

    def read_tree_file(self):
        return Phylo.read(self.tree_file_path, self.tree_format)

    def read_tree1_file(self):
        return Phylo.read(self.tree1_file_path, self.tree_format)

    def write_tree_file(self, tree, output_file_path):
        return Phylo.write(tree, output_file_path, self.tree_format)
