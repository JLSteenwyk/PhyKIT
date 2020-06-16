import getopt
import logging
import os.path
import sys

from scipy.stats import zscore
from scipy.stats.stats import pearsonr

from Bio import Phylo

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
        tree_zero_tips = self.tip_names(tree_zero)
        tree_one_tips = self.tip_names(tree_one)
        tree_ref_tips = self.tip_names(tree_ref)

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

        # obtain corrected branch lengths where branch lengths
        # are corrected by the species tree branch length
        tree_zero_corr_branch_lengths = self.correct_branch_lengths(tree_zero, tree_ref)
        tree_one_corr_branch_lengths = self.correct_branch_lengths(tree_one, tree_ref)

        # remove corrected BLs greater than 5
        outlier_indices = []
        for idx in range(0, len(tree_zero_corr_branch_lengths)): 
            if tree_zero_corr_branch_lengths[idx] > 5 or tree_zero_corr_branch_lengths[idx] < -5: 
                outlier_indices.append(idx) 
        for idx in range(0, len(tree_one_corr_branch_lengths)): 
            if tree_one_corr_branch_lengths[idx] > 5 or tree_one_corr_branch_lengths[idx] < -5:
                if idx not in outlier_indices: 
                    outlier_indices.append(idx)  
        tree_zero_corr_branch_lengths = [i for j, i in enumerate(tree_zero_corr_branch_lengths) if j not in outlier_indices]
        tree_one_corr_branch_lengths = [i for j, i in enumerate(tree_one_corr_branch_lengths) if j not in outlier_indices]

        corrected_branch_lengths = []
        temp_branch_lengths = []
        temp_branch_lengths.append(tree_zero_corr_branch_lengths)
        temp_branch_lengths.append(tree_one_corr_branch_lengths)
        corrected_branch_lengths.append(temp_branch_lengths) 

        # standardize values for final correction
        tree_zero_corr_branch_lengths = zscore(tree_zero_corr_branch_lengths)
        tree_one_corr_branch_lengths = zscore(tree_one_corr_branch_lengths)

        # Calculate correlation and append to results array
        # also keep a list of p values
        corr = (list(pearsonr(tree_zero_corr_branch_lengths, tree_one_corr_branch_lengths)))
        
        if self.verbose:
            for val_zero, val_one in zip(tree_zero_corr_branch_lengths, tree_one_corr_branch_lengths):
                print(val_zero, val_one)
        else:
            print(corr[0], corr[1])

    def process_args(self, args):
        return dict(
            tree_file_path=args.tree_zero,
            tree1_file_path=args.tree_one,
            reference=args.reference,
            verbose=args.verbose
        )

    def tip_names(
        self, tree
        ):
        """
        get tip names from a tree
        """

        tips = []
        for tip in tree.get_terminals():
            tips.append(tip.name)
        return(tips)

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
    
    def prune_tips(
        self, tree, tips
        ):
        """
        prune tips from trees
        """

        for tip in tips:
            tree.prune(tip)
        return(tree)

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