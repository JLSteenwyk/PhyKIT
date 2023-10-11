import copy
import sys

from scipy.stats import (pearsonr, zscore)

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
        tree_zero_tips_to_prune = list(set(tree_zero_tips) - set(shared_tree_tips))
        tree_one_tips_to_prune = list(set(tree_one_tips) - set(shared_tree_tips))
        tree_ref_tips_to_prune = list(set(tree_ref_tips) - set(shared_tree_tips))

        # get a set of pruned trees
        tree_zero = self.prune_tips(tree_zero, tree_zero_tips_to_prune)
        tree_one = self.prune_tips(tree_one, tree_one_tips_to_prune)
        tree_ref = self.prune_tips(tree_ref, tree_ref_tips_to_prune)

        # obtain corrected branch lengths where branch lengths
        # are corrected by the species tree branch length
        tree_zero_corr_branch_lengths, tree_one_corr_branch_lengths, tip_names = self.correct_branch_lengths(tree_zero, tree_one, tree_ref)

        # remove corrected BLs greater than 5
        outlier_indices = []
        outlier_indices = self.get_indices_of_outlier_branch_lengths(tree_zero_corr_branch_lengths, outlier_indices)
        outlier_indices = self.get_indices_of_outlier_branch_lengths(tree_one_corr_branch_lengths, outlier_indices)

        tree_zero_corr_branch_lengths = self.remove_outliers_based_on_indices(tree_zero_corr_branch_lengths, outlier_indices)
        tree_one_corr_branch_lengths = self.remove_outliers_based_on_indices(tree_one_corr_branch_lengths, outlier_indices)
        tip_names = self.remove_outliers_based_on_indices(tip_names, outlier_indices)

        # standardize values for final correction
        tree_zero_corr_branch_lengths = zscore(tree_zero_corr_branch_lengths)
        tree_one_corr_branch_lengths = zscore(tree_one_corr_branch_lengths)

        # Calculate correlation and append to results array
        # also keep a list of p values
        corr = (list(pearsonr(tree_zero_corr_branch_lengths, tree_one_corr_branch_lengths)))

        try:
            if self.verbose:
                for val_zero, val_one, tip_name in zip(tree_zero_corr_branch_lengths, tree_one_corr_branch_lengths, tip_names):
                    print(f"{round(val_zero, 4)}\t{round(val_one, 4)}\t{';'.join(tip_name)}")
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

    def get_indices_of_outlier_branch_lengths(
        self, corr_branch_lengths, outlier_indices
    ):
        """
        create index for branch lengths that 
        have an absolute value greater than 5
        """
        for idx in range(0, len(corr_branch_lengths)): 
            try:
                if corr_branch_lengths[idx] > 5 or corr_branch_lengths[idx] < -5:
                    if idx not in outlier_indices: 
                        outlier_indices.append(idx)
            except TypeError:
                outlier_indices.append(idx)

        return outlier_indices

    def remove_outliers_based_on_indices(self, corr_branch_lengths, outlier_indices):
        """
        remove value if the value is an outlier according
        to the outlier indices list
        """
        corr_branch_lengths = [i for j, i in enumerate(corr_branch_lengths) if j not in outlier_indices]
        return corr_branch_lengths

    def prune_tips(
        self, tree, tips
    ):
        """
        prune tips from trees
        """

        for tip in tips:
            tree.prune(tip)

        return tree

    def correct_branch_lengths(
        self,
        t0,
        t1,
        sp
    ):
        """
        obtain a list of corrected branch lengths
        """

        l0 = []
        l1 = []
        tip_names = []
        # terminal corrected branch lengths
        for i in sp.get_terminals():
            newtree = copy.deepcopy(t0)
            newtree1 = copy.deepcopy(t1)
            sp_tips = self.get_tip_names_from_tree(i)
            tip_names.append(sp_tips)
            newtree = newtree.common_ancestor(i.name)
            newtree1 = newtree1.common_ancestor(i.name)
            l0.append(round(newtree.branch_length / i.branch_length, 6))
            l1.append(round(newtree1.branch_length / i.branch_length, 6))

        # nonterminal corrected branch lengths
        for i in sp.get_nonterminals():
            newtree = copy.deepcopy(t0)
            newtree1 = copy.deepcopy(t1)
            sp_tips = self.get_tip_names_from_tree(i)            
            newtree = newtree.common_ancestor(sp_tips)
            newtree1 = newtree1.common_ancestor(sp_tips)
            try:
                l0.append(round(newtree.branch_length / i.branch_length, 6))
                tip_names.append(sp_tips)
            except TypeError:
                continue
            try:
                l1.append(round(newtree1.branch_length / i.branch_length, 6))
            except TypeError:
                continue

        return (l0, l1, tip_names)
