import copy
import statistics as stat
import sys

import numpy as np

from .base import Tree


class HiddenParalogyCheck(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        master_tree = self.read_tree_file()
        clades = self.read_clades_file(self.clade)

        # initialize results array
        res_arr = []

        for clade in clades:
            tree = copy.deepcopy(master_tree)
            clan = [{"name": taxon_name} for taxon_name in clade]

            # get tip names
            tree_tips = self.get_tip_names_from_tree(tree)
            clade_of_interest = list(set(clade).intersection(tree_tips))

            # check for sufficient representation
            if len(clade_of_interest) <= 1:
                res_arr.append(["insufficient_taxon_representation"])
                continue

            # determine shared tips
            shared_tree_tips = self.shared_tips(clade_of_interest, tree_tips)
            # get tips that differ
            diff_tips = list(set(tree_tips) - set(shared_tree_tips))
            # root tree with different tips
            tree.root_with_outgroup(diff_tips)
            # get common ancestor of clan of interest
            tree = tree.common_ancestor(shared_tree_tips)
            # get common ancestor tree tips
            common_ancestor_tips = self.get_tip_names_from_tree(tree)
            diff_tips_between_clade_and_curr_tree = list(
                set(clade_of_interest) ^ set(common_ancestor_tips)
            )

            res_arr = self.populate_res_arr(
                shared_tree_tips, diff_tips_between_clade_and_curr_tree, res_arr
            )

        self.print_results(res_arr)

    def process_args(self, args):
        tree_file_path = args.tree

        return dict(
            tree_file_path=tree_file_path,
            clade=args.clade,
        )

    def read_clades_file(self, clades):
        try:
            clades = [[s for s in l.split()] for l in open(self.clade).readlines()]
        except FileNotFoundError:
            print("Clade file is not found. Please check pathing.")
            sys.exit()

        return clades

    def populate_res_arr(
        self, shared_tree_tips, diff_tips_between_clade_and_curr_tree, res_arr
    ):
        temp = []

        if len(diff_tips_between_clade_and_curr_tree) == 0:
            temp.append("monophyletic")
        else:
            temp.append("not_monophyletic")
        temp.append(diff_tips_between_clade_and_curr_tree)
        res_arr.append(temp)

        return res_arr

    def print_results(self, res_arr):
        for res in res_arr:
            print(
                f"{res[0]}"
            )
