import math
from typing import Dict, Tuple

from Bio.Phylo import Newick

from .base import Tree
from ...helpers.json_output import print_json


class KuhnerFelsensteinDistance(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            tree1_file_path=parsed["tree1_file_path"],
        )
        self.json_output = parsed["json_output"]

    def run(self):
        tree_zero = self.read_tree_file()
        tree_one = self.read_tree1_file()

        # Get shared tree tip names
        tree_zero_tips = set(self.get_tip_names_from_tree(tree_zero))
        tree_one_tips = set(self.get_tip_names_from_tree(tree_one))
        shared_tree_tips = tree_zero_tips & tree_one_tips

        if not shared_tree_tips:
            print("Trees share no common taxa.")
            raise SystemExit(2)

        # Prune to common set
        tree_zero_tips_to_prune = list(tree_zero_tips - shared_tree_tips)
        tree_one_tips_to_prune = list(tree_one_tips - shared_tree_tips)

        if tree_zero_tips_to_prune:
            tree_zero = self.prune_tree_using_taxa_list(tree_zero, tree_zero_tips_to_prune)
        if tree_one_tips_to_prune:
            tree_one = self.prune_tree_using_taxa_list(tree_one, tree_one_tips_to_prune)

        # Root on same taxon
        tip_for_rooting = tree_zero.get_terminals()[0].name
        tree_zero.root_with_outgroup(tip_for_rooting)
        tree_one.root_with_outgroup(tip_for_rooting)

        plain_kf, normalized_kf = self.calculate_kf_distance(tree_zero, tree_one)

        if self.json_output:
            print_json(
                dict(
                    plain_kf=round(plain_kf, 4),
                    normalized_kf=round(normalized_kf, 4),
                )
            )
            return

        print(f"{round(plain_kf, 4)}\t{round(normalized_kf, 4)}")

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            tree_file_path=args.tree_zero,
            tree1_file_path=args.tree_one,
            json_output=getattr(args, "json", False),
        )

    def _get_splits_with_lengths(self, tree) -> Dict[frozenset, float]:
        """Return dict mapping frozenset(tip_names) -> branch_length.

        Includes both internal and terminal branches.
        None branch lengths are treated as 0.0.
        The root node itself is excluded (it has no parent branch).
        """
        splits = {}
        for clade in tree.find_clades(order="postorder"):
            if clade == tree.root:
                continue
            tips = frozenset(t.name for t in clade.get_terminals())
            bl = clade.branch_length if clade.branch_length is not None else 0.0
            splits[tips] = bl
        return splits

    def calculate_kf_distance(
        self, tree_zero: Newick.Tree, tree_one: Newick.Tree
    ) -> Tuple[float, float]:
        splits_zero = self._get_splits_with_lengths(tree_zero)
        splits_one = self._get_splits_with_lengths(tree_one)

        all_splits = set(splits_zero.keys()) | set(splits_one.keys())

        kf_squared = 0.0
        for split in all_splits:
            b0 = splits_zero.get(split, 0.0)
            b1 = splits_one.get(split, 0.0)
            kf_squared += (b0 - b1) ** 2

        plain_kf = math.sqrt(kf_squared)

        total_bl = sum(splits_zero.values()) + sum(splits_one.values())
        normalized_kf = plain_kf / total_bl if total_bl > 0 else 0.0

        return plain_kf, normalized_kf
