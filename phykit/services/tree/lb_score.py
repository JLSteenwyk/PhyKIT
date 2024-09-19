import sys
import itertools
from typing import Dict, List, Tuple

from Bio.Phylo import Newick

from .base import Tree

from ...helpers.stats_summary import (
    calculate_summary_statistics_from_arr,
    print_summary_statistics,
)


class LBScore(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        tree = self.read_tree_file()
        tips, LBis = self.calculate_lb_score(tree)
        if self.verbose:
            try:
                for tip, LBi in zip(tips, LBis):
                    print(f"{tip}\t{round(LBi, 4)}")
            except BrokenPipeError:
                pass
        else:
            stats = calculate_summary_statistics_from_arr(LBis)
            print_summary_statistics(stats)

    def process_args(self, args) -> Dict[str, str]:
        return dict(tree_file_path=args.tree, verbose=args.verbose)

    def calculate_average_distance_between_tips(
        self,
        tips: List[str],
        tree: Newick.Tree,
    ) -> float:
        total_dist = sum(
            tree.distance(tip1, tip2)
            for tip1, tip2 in itertools.combinations(tips, 2)
        )

        num_combos = len(tips) * (len(tips) - 1) // 2

        return total_dist / num_combos if num_combos else 0

    def calculate_average_distance_of_taxon_to_other_taxa(
        self,
        tips: List[str],
        tree: Newick.Tree,
    ) -> List[float]:
        avg_PDis = []
        for tip in tips:
            tips_minus_i = list(set(tips) - set(tip))
            PDi = []
            for tip_minus in tips_minus_i:
                PDi.append(tree.distance(tip, tip_minus))
            PDi = sum(PDi) / len(PDi)
            avg_PDis.append(PDi)

        return avg_PDis

    def calculate_lb_score_per_taxa(
        self,
        avg_PDis: List[float],
        avg_dist: float
    ) -> List[float]:
        LBis = []
        for PDi in avg_PDis:
            try:
                LBis.append((((PDi / avg_dist) - 1) * 100))
            except ZeroDivisionError:
                try:
                    print("Invalid tree. Tree should contain branch lengths")
                    sys.exit()
                except BrokenPipeError:
                    pass

        return LBis

    def calculate_lb_score(
        self,
        tree: Newick.Tree
    ) -> Tuple[List[str], List[float]]:
        tips = self.get_tip_names_from_tree(tree)

        avg_dist = self.calculate_average_distance_between_tips(tips, tree)

        avg_PDis = \
            self.calculate_average_distance_of_taxon_to_other_taxa(tips, tree)

        LBis = self.calculate_lb_score_per_taxa(avg_PDis, avg_dist)

        return tips, LBis
