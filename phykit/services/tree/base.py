from Bio import Phylo

from ..base import BaseService

class Tree(BaseService):
    def __init__(
        self,
        *args,
        tree_file_path,
        alignment_file_path=None,
        tree1_file_path=None,
        outgroup_taxa_file_path=None,
        output_file_path=None,
        factor=None,
        verbose=None
    ):
        self.tree_file_path = tree_file_path
        self.tree1_file_path = tree1_file_path
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

    def calculate_treeness(self, tree=None, print_value=False):
        if not tree:
            tree = self.read_tree_file()

        inter_len = float(0.0)
        # determine internal branch lengths
        for interal in tree.get_nonterminals():
            # only include if a branch length value is present
            if interal.branch_length != None:
                inter_len += interal.branch_length
        # determine total branch length
        total_len = tree.total_branch_length()

        try:
            treeness = float(inter_len / total_len)
            if print_value:
                print(f"{treeness}")
            return treeness
        except ZeroDivisionError:
            print("Invalid tree. Tree should contain branch lengths")
            return None
