from typing import Dict, List

from Bio.Phylo import Newick

from .base import Tree


class RobinsonFouldsDistance(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree_zero = self.read_tree_file()
        tree_one = self.read_tree1_file()

        # get shared tree tip names
        tree_zero_tips = self.get_tip_names_from_tree(tree_zero)
        tree_one_tips = self.get_tip_names_from_tree(tree_one)
        shared_tree_tips = self.shared_tips(tree_zero_tips, tree_one_tips)

        # prune to common set
        tree_zero_tips_to_prune = list(set(tree_zero_tips) - set(shared_tree_tips))
        tree_one_tips_to_prune = list(set(tree_one_tips) - set(shared_tree_tips))
        tree_zero = self.prune_tree_using_taxa_list(tree_zero, tree_zero_tips_to_prune)
        tree_one = self.prune_tree_using_taxa_list(tree_one, tree_one_tips_to_prune)

        for term in tree_zero.get_terminals():
            tip_for_rooting = term.name
            break
        tree_zero.root_with_outgroup(tip_for_rooting)
        tree_one.root_with_outgroup(tip_for_rooting)

        plain_rf, normalized_rf = self.calculate_robinson_foulds_distance(
            tree_zero, tree_one
        )

        print(f"{plain_rf}\t{round(normalized_rf, 4)}")

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            tree_file_path=args.tree_zero,
            tree1_file_path=args.tree_one,
        )

    def calculate_robinson_foulds_distance(self, tree_zero, tree_one):
        plain_rf = 0
        plain_rf = self.compare_trees(plain_rf, tree_zero, tree_one)
        plain_rf = self.compare_trees(plain_rf, tree_one, tree_zero)

        tip_count = tree_zero.count_terminals()

        normalized_rf = plain_rf / (2 * (tip_count - 3))

        return plain_rf, normalized_rf

    def compare_trees(
        self,
        plain_rf: int,
        tree_zero: Newick.Tree,
        tree_one: Newick.Tree
    ) -> int:
        # loop through tree_zero and find similar clade in tree_one
        for clade_zero in tree_zero.get_nonterminals()[1:]:
            # initialize and populate a list of tip names in tree_zero
            tip_names_zero = self.get_tip_names_from_tree(clade_zero)
            # get common ancestor of tree_zero tip names in tree_one
            clade_one = tree_one.common_ancestor(tip_names_zero)
            # initialize and populate a list of tip names in tree_one
            tip_names_one = self.get_tip_names_from_tree(clade_one)
            # compare the list of tip names
            plain_rf = self.determine_if_clade_differs(
                plain_rf, tip_names_zero, tip_names_one
            )

        return plain_rf

    def determine_if_clade_differs(
        self,
        plain_rf: int,
        tip_names_zero: List[str],
        tip_names_one: List[str],
    ) -> int:
        if set(tip_names_zero) != set(tip_names_one):
            plain_rf += 1

        return plain_rf
