import getopt
import logging
import os.path
import sys

from scipy.stats import zscore
from scipy.stats.stats import pearsonr

from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin

from .base import Tree

class CovaryingEvolutionaryRates(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree_zero = self.read_tree_file()
        tree_one = self.read_tree1_file()
        tree_ref = self.read_reference_tree_file()

        ## - Calculate correlation between two gene trees 
        ## and save results to an array, corrArr.
        ## - Branch lengths will also be part of output

        # get tree tip names
        tree_zero_tips = self.get_tip_names_from_tree(tree_zero)
        tree_one_tips = self.get_tip_names_from_tree(tree_one)
        tree_ref_tips = self.get_tip_names_from_tree(tree_ref)

        # get shared tips between the two trees
        shared_tree_tips = self.shared_tips(tree_zero_tips, tree_one_tips)

        # find differences between tree tips and shared tips
        # to determine what tips to prune
        tree_zero_tips_to_prune   = list(set(tree_zero_tips) - set(shared_tree_tips))
        tree_one_tips_to_prune    = list(set(tree_one_tips) - set(shared_tree_tips))
        tree_ref_tips_to_prune    = list(set(tree_ref_tips) - set(shared_tree_tips))

        # get a set of pruned trees
        tree_zero = self.prune_tips(tree_zero, tree_zero_tips_to_prune)
        tree_one = self.prune_tips(tree_one, tree_one_tips_to_prune)
        tree_ref = self.prune_tips(tree_ref, tree_ref_tips_to_prune)

        # check that the input trees have the same topology
        self.determine_if_trees_differ(tree_zero, tree_one, tree_ref)

        # obtain corrected branch lengths where branch lengths
        # are corrected by the species tree branch length
        tree_zero_corr_branch_lengths = self.correct_branch_lengths(tree_zero, tree_ref)
        tree_one_corr_branch_lengths = self.correct_branch_lengths(tree_one, tree_ref)

        # remove corrected BLs greater than 5
        outlier_indices = []
        outlier_indices = self.get_indices_of_outlier_branch_lengths(tree_zero_corr_branch_lengths, outlier_indices)
        outlier_indices = self.get_indices_of_outlier_branch_lengths(tree_one_corr_branch_lengths, outlier_indices)

        tree_zero_corr_branch_lengths = self.remove_outliers_based_on_indices(tree_zero_corr_branch_lengths, outlier_indices)
        tree_one_corr_branch_lengths = self.remove_outliers_based_on_indices(tree_one_corr_branch_lengths, outlier_indices)

        # standardize values for final correction
        tree_zero_corr_branch_lengths = zscore(tree_zero_corr_branch_lengths)
        tree_one_corr_branch_lengths = zscore(tree_one_corr_branch_lengths)

        # Calculate correlation and append to results array
        # also keep a list of p values
        corr = (list(pearsonr(tree_zero_corr_branch_lengths, tree_one_corr_branch_lengths)))
        
        try:
            if self.verbose:
                for val_zero, val_one in zip(tree_zero_corr_branch_lengths, tree_one_corr_branch_lengths):
                    print(f"{round(val_zero, 4)}\t{round(val_one, 4)}")
            else:
                print(f"{round(corr[0], 4)}\t{round(corr[1], 6)}")
        except BrokenPipeError:
            pass

    def process_args(self, args):
        return dict(
            tree_file_path=args.tree_zero,
            tree1_file_path=args.tree_one,
            reference=args.reference,
            verbose=args.verbose
        )

    def get_indices_of_outlier_branch_lengths(self, corr_branch_lengths, outlier_indices):
        """
        create index for branch lengths that 
        have an absolute value greater than 5
        """
        for idx in range(0, len(corr_branch_lengths)): 
            if corr_branch_lengths[idx] > 5 or corr_branch_lengths[idx] < -5:
                if idx not in outlier_indices: 
                    outlier_indices.append(idx)
        
        return outlier_indices

    def remove_outliers_based_on_indices(self, corr_branch_lengths, outlier_indices):
        """
        remove value if the value is an outlier according
        to the outlier indices list
        """
        corr_branch_lengths = [i for j, i in enumerate(corr_branch_lengths) if j not in outlier_indices]
        return corr_branch_lengths

    # def shared_tips(
    #     self, 
    #     a,
    #     b
    #     ):
    #     """
    #     Determines what tips are shared between two trees
    #     -------------------------------------------------
    #     argv: a
    #         list of tips from one tree
    #     argv: b
    #         list of tips from a second tree
    #     """ 

    #     a_set = set(a) 
    #     b_set = set(b) 
        
    #     # check length  
    #     if len(a_set.intersection(b_set)) > 0: 
    #         return(list(a_set.intersection(b_set)))   
    #     else: 
    #         print("no common tips") 
    #         sys.exit()
 
    
    def prune_tips(
        self, tree, tips
        ):
        """
        prune tips from trees
        """

        for tip in tips:
            tree.prune(tip)
        
        return tree

    def determine_if_trees_differ(self, tree_zero, tree_one, tree_ref):
        """
        determine if the trees differ from one another
        """
        differences = 0
        differences = self.compare_trees(differences, tree_zero, tree_one)
        differences = self.compare_trees(differences, tree_zero, tree_ref)
        differences = self.compare_trees(differences, tree_one, tree_ref)
        
        if differences > 0:
            print("Input trees differ in topology.")
            print("Please ensure input phylogenies all have the same topology.")
            sys.exit()

    def compare_trees(
        self,
        plain_rf: int,
        tree_zero: Tree,
        tree_one: Tree
    ) -> int:
        # loop through tree_zero and find similar clade in tree_one
        for clade_zero in tree_zero.get_nonterminals():
            # initialize and populate a list of tip names in tree_zero
            tip_names_zero = self.get_tip_names_from_tree(clade_zero)
            # get common ancestor of tree_zero tip names in tree_one
            clade_one = tree_one.common_ancestor(tip_names_zero)
            # initialize and populate a list of tip names in tree_one
            tip_names_one = self.get_tip_names_from_tree(clade_one)

            # compare the list of tip names
            plain_rf = self.determine_if_clade_differs(
                plain_rf,
                tip_names_zero,
                tip_names_one
            )

        return plain_rf
            
    def determine_if_clade_differs(
        self,
        plain_rf: int,
        tip_names_zero: list,
        tip_names_one: list
        ) -> int:
        """
        if clade differs, add 1 to plain_rf value
        """

        if set(tip_names_zero) != set(tip_names_one):
            plain_rf +=1
        
        return plain_rf

    def correct_branch_lengths(
        self, 
        t,
        sp
        ):
        """
        obtain a list of corrected branch lengths
        """
        
        l = []
        for i, s in zip(t.get_terminals(), sp.get_terminals()):
            l.append(i.branch_length/s.branch_length)
        return(l)