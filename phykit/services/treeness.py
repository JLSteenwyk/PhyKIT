from Bio import Phylo
from phykit.services.base import BaseService


class Treeness(BaseService):
    def __init__(self, args):
        self.process_args(args)

    def run(self):
        tree = self.read_file()
        treeness = self.calculate_treeness(tree)
        if treeness:
            print(f"Treeness score: {treeness}")

    def process_args(self, args):
        self.tree_file_path = args.tree

    def read_file(self):
        return Phylo.read(self.tree_file_path, "newick")

    def calculate_treeness(self, tree):
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
            return treeness
        except ZeroDivisionError:
            print("Invalid tree. Tree should contain branch lengths")
            return None
