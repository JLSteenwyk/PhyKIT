import copy
import sys

from typing import Dict, List, Union

from .base import Tree


class HiddenParalogyCheck(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        master_tree = self.read_tree_file()
        clades = self.read_clades_file(self.clade)

        res_arr = []

        master_tree_tips = set(self.get_tip_names_from_tree(master_tree))

        for clade in clades:
            master_tree = copy.deepcopy(master_tree)
            clade_of_interest = set(clade).intersection(master_tree_tips)

            if len(clade_of_interest) <= 1:
                res_arr.append(["insufficient_taxon_representation"])
                continue

            diff_tips = master_tree_tips - clade_of_interest
            master_tree.root_with_outgroup(list(diff_tips))

            subtree = master_tree.common_ancestor(clade_of_interest)
            common_ancestor_tips = set(self.get_tip_names_from_tree(subtree))

            diff_tips_between_clade_and_curr_tree = \
                clade_of_interest.symmetric_difference(common_ancestor_tips)

            res_arr.append(
                [
                    "monophyletic" if not diff_tips_between_clade_and_curr_tree else "not_monophyletic",
                    list(diff_tips_between_clade_and_curr_tree),
                ]
            )

        self.print_results(res_arr)

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            tree_file_path=args.tree,
            clade=args.clade,
        )

    def read_clades_file(self, clades: str) -> List[List[str]]:
        try:
            with open(clades, 'r') as file:
                return [line.split() for line in file.readlines()]
        except FileNotFoundError:
            print("Clade file not found. Please check the path.")
            sys.exit()

    def print_results(self, res_arr: List[List[Union[List, str]]]) -> None:
        for res in res_arr:
            print(res[0])
