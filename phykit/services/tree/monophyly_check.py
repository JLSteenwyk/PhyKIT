import sys

from Bio import Phylo

from .base import Tree

from ...helpers.stats_summary import calculate_summary_statistics_from_arr
from ...helpers.files import read_single_column_file_to_list


class MonophylyCheck(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        taxa = read_single_column_file_to_list(self.list_of_taxa)

        # initialize results array
        res_arr = []

        clan = [{"name": taxon_name} for taxon_name in taxa]

        # get tip names
        tree_tips = self.get_tip_names_from_tree(tree)
        taxa_of_interest = list(set(taxa).intersection(tree_tips))

        # check for sufficient representation
        if len(taxa_of_interest)<=1:
            res_arr.append(["insufficient_taxon_representation"])
            sys.exit()

        # determine shared tips
        shared_tree_tips = self.shared_tips(taxa_of_interest, tree_tips)
        # get tips that differ
        diff_tips = list(set(tree_tips) - set(shared_tree_tips))
        # root tree with different tips
        tree.root_with_outgroup(diff_tips)
        # get common ancestor of clan of interest
        tree = tree.common_ancestor(shared_tree_tips)
        # get common ancestor tree tips
        common_ancestor_tips = self.get_tip_names_from_tree(tree)        
        diff_tips_between_clade_and_curr_tree = list(set(taxa_of_interest) ^ set(common_ancestor_tips))

        stats = self.get_bootstrap_statistics(tree)

        res_arr = self.populate_res_arr(
            diff_tips_between_clade_and_curr_tree,
            stats,
            res_arr
        )

        self.print_results(res_arr)

    def process_args(self, args):
        return dict(
            tree_file_path=args.tree,
            list_of_taxa=args.list_of_taxa,
        )
    
    def get_bootstrap_statistics(self, tree):
        # get bootstrap support values
        bs_vals = []
        # populate bs_vals with bootstrap values
        for terminal in tree.get_nonterminals():
        # only include if a bootstrap value is present
            if terminal.confidence != None:
                bs_vals.append(terminal.confidence)
        stats = calculate_summary_statistics_from_arr(bs_vals)

        return stats

    def populate_res_arr(
        self,
        diff_tips_between_clade_and_curr_tree,
        stats,
        res_arr
    ):

        temp = []

        if len(diff_tips_between_clade_and_curr_tree) == 0:
            temp.append("monophyletic")
        else:
            temp.append("not_monophyletic")
        temp.append(stats["mean"])
        temp.append(stats["maximum"])
        temp.append(stats["minimum"])
        temp.append(stats["standard_deviation"])
        temp.append(diff_tips_between_clade_and_curr_tree)
        res_arr.append(temp)

        return res_arr

    def print_results(self, res_arr):
        for res in res_arr:
            try:
                if len(res[5]) != 0:
                    res[5].sort()
                    print(f"{res[0]}\t{round(res[1], 4)}\t{round(res[2], 4)}\t{round(res[3], 4)}\t{round(res[4], 4)}\t{';'.join(res[5])}")
                else:
                    print(f"{res[0]}\t{round(res[1], 4)}\t{round(res[2], 4)}\t{round(res[3], 4)}\t{round(res[4], 4)}")
            except IndexError:
                print(f"{res[0]}")