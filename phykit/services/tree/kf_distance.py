from __future__ import annotations

import math

from .base import Tree


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class KuhnerFelsensteinDistance(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            tree1_file_path=parsed["tree1_file_path"],
        )
        self.json_output = parsed["json_output"]

    def run(self):
        if self.tree_file_path == self.tree1_file_path:
            self._output_result(0.0, 0.0)
            return

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
        tip_for_rooting = self.get_first_tip_name_from_tree(tree_zero)
        tree_zero.root_with_outgroup(tip_for_rooting)
        tree_one.root_with_outgroup(tip_for_rooting)

        plain_kf, normalized_kf = self.calculate_kf_distance(tree_zero, tree_one)

        self._output_result(plain_kf, normalized_kf)

    def _output_result(self, plain_kf, normalized_kf):
        if self.json_output:
            print_json(
                dict(
                    plain_kf=round(plain_kf, 4),
                    normalized_kf=round(normalized_kf, 4),
                )
            )
            return

        print(f"{round(plain_kf, 4)}\t{round(normalized_kf, 4)}")

    def process_args(self, args) -> dict[str, str]:
        return dict(
            tree_file_path=args.tree_zero,
            tree1_file_path=args.tree_one,
            json_output=getattr(args, "json", False),
        )

    def _get_splits_with_lengths(self, tree) -> dict[frozenset, float]:
        """Return dict mapping frozenset(tip_names) -> branch_length.

        Includes both internal and terminal branches.
        None branch lengths are treated as 0.0.
        The root node itself is excluded (it has no parent branch).
        """
        postorder = self._postorder_clades_direct(tree)
        if postorder is None:
            postorder = tree.find_clades(order="postorder")

        clade_tips = {}
        splits = {}
        root = tree.root
        for clade in postorder:
            children = getattr(clade, "clades", None)
            if children:
                child_count = len(children)
                if child_count == 2:
                    tip_set = (
                        clade_tips[id(children[0])]
                        | clade_tips[id(children[1])]
                    )
                elif child_count == 1:
                    tip_set = clade_tips[id(children[0])]
                else:
                    tip_set = frozenset().union(
                        *(clade_tips[id(child)] for child in children)
                    )
            else:
                tip_set = frozenset([clade.name])

            clade_tips[id(clade)] = tip_set
            if clade is not root:
                bl = (
                    clade.branch_length
                    if clade.branch_length is not None
                    else 0.0
                )
                splits[tip_set] = bl
        return splits

    @staticmethod
    def _postorder_clades_direct(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        clades = []
        stack = [root]
        append_clade = clades.append
        pop = stack.pop
        extend = stack.extend
        try:
            while stack:
                clade = pop()
                append_clade(clade)
                extend(clade.clades)
        except AttributeError:
            return None
        clades.reverse()
        return clades

    def calculate_kf_distance(
        self, tree_zero: Newick.Tree, tree_one: Newick.Tree
    ) -> tuple[float, float]:
        if tree_zero is tree_one:
            return 0.0, 0.0

        splits_zero = self._get_splits_with_lengths(tree_zero)
        splits_one = self._get_splits_with_lengths(tree_one)

        kf_squared = 0.0
        total_bl = 0.0
        remaining_splits_one = splits_one.copy()
        for split, b0 in splits_zero.items():
            b1 = remaining_splits_one.pop(split, 0.0)
            diff = b0 - b1
            kf_squared += diff * diff
            total_bl += b0 + b1

        for b1 in remaining_splits_one.values():
            kf_squared += b1 * b1
            total_bl += b1

        plain_kf = math.sqrt(kf_squared)
        normalized_kf = plain_kf / total_bl if total_bl > 0 else 0.0

        return plain_kf, normalized_kf
