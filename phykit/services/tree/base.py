import sys
from typing import Tuple

from Bio import Phylo

from ..base import BaseService

class Tree(BaseService):
    def __init__(
        self,
        *args,
        tree_file_path=None,
        idmap=None,
        alignment_file_path=None,
        tree1_file_path=None,
        outgroup_taxa_file_path=None,
        output_file_path=None,
        factor=None,
        remove=None,
        verbose=None,
        reference=None,
        list_of_taxa=None,
        trees=None,
        groups=None,
        support=None,
        tip_1=None,
        tip_2=None,
        clade=None,
        keep=None,
    ):
        self.tree_file_path = tree_file_path
        self.tree1_file_path = tree1_file_path
        self.alignment_file_path = alignment_file_path
        self.output_file_path = output_file_path
        self.outgroup_taxa_file_path = outgroup_taxa_file_path
        self.tree_format = "newick"
        self.verbose = verbose
        self.factor = factor
        self.remove = remove
        self.idmap = idmap
        self.reference = reference
        self.list_of_taxa = list_of_taxa
        self.trees = trees
        self.groups = groups
        self.support = support
        self.tip_1 = tip_1
        self.tip_2 = tip_2
        self.clade = clade
        self.keep = keep

    def read_tree_file(self):
        try:
            return Phylo.read(self.tree_file_path, self.tree_format)
        except FileNotFoundError:
            print(f"{self.tree_file_path} corresponds to no such file or directory.")
            print("Please check filename and pathing")
            sys.exit()


    def read_tree1_file(self):
        try:
            return Phylo.read(self.tree1_file_path, self.tree_format)
        except FileNotFoundError:
            print(f"{self.tree1_file_path} corresponds to no such file or directory.")
            print("Please check filename and pathing")
            sys.exit()


    def read_reference_tree_file(self):
        try:
            return Phylo.read(self.reference, self.tree_format)
        except FileNotFoundError:
            print(f"{self.reference} corresponds to no such file or directory.")
            print("Please check filename and pathing")
            sys.exit()


    def write_tree_file(self, tree, output_file_path):
        return Phylo.write(tree, output_file_path, self.tree_format)

    def get_tip_names_from_tree(
        self, tree
        ) -> list:
        """
        get tip names from a tree
        """
        tips = []
        for tip in tree.get_terminals():
            tips.append(tip.name)
        
        return tips

    def shared_tips(
        self, 
        a,
        b
        ):
        """
        Determines what tips are shared between two trees
        -------------------------------------------------
        argv: a
            list of tips from one tree
        argv: b
            list of tips from a second tree
        """ 

        a_set = set(a) 
        b_set = set(b) 
        
        # check length  
        if len(a_set.intersection(b_set)) > 0: 
            return(list(a_set.intersection(b_set)))   
        else: 
            print("no common tips") 
            sys.exit()

    def prune_tree_using_taxa_list(
        self,
        tree,
        taxa_to_prune: list
    ):
        """
        prune taxa from tree
        """
        for taxon in taxa_to_prune:
            tree.prune(taxon)
        
        return tree

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
            try:
                if print_value:
                    print(f"{treeness}")
                return treeness
            except BrokenPipeError:
                pass
        except ZeroDivisionError:
            try:
                print("Invalid tree. Tree should contain branch lengths")
                return None
            except BrokenPipeError:
                pass
