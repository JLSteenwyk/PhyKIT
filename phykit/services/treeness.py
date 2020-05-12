from phykit.services.base import BaseService

from Bio import Phylo

class Treeness(BaseService):
    def __init__(self, args):
        self.process_args(args)

    def run(self):
        tree = self.read_file()
        treeness = self.calculate_treeness(tree)

    def process_args(self, args):
        self.tree_file_path = args.tree
    
    def read_file(self):
        return Phylo.read(self.tree_file_path, 'newick')

    def calculate_treeness(self, tree):
        # initialize variables for terminal branch length
        interLen = float(0.0)
        # determine internal branch lengths
        for interal in tree.get_nonterminals():
            # only include if a branch length value is present
            if interal.branch_length != None:
                interLen += interal.branch_length

        # initialize variable for total branch length
        totalLen = float(0.0)
        # determine total branch length
        totalLen = tree.total_branch_length()

        return float(interLen/totalLen)

