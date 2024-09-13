from typing import Dict, List, Tuple
import itertools

from Bio.Phylo import Newick

from .base import Tree

from ...helpers.stats_summary import (
    calculate_summary_statistics_from_arr,
    print_summary_statistics,
)


class PatristicDistances(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        patristic_distances, combos, stats = \
            self.calculate_patristic_distances(tree)

        if self.verbose:
            try:
                for combo, patristic_distance in zip(combos, patristic_distances):
                    print(f"{combo[0]}\t{combo[1]}\t{round(patristic_distance, 4)}")
            except BrokenPipeError:
                pass
        else:
            print_summary_statistics(stats)

    def process_args(self, args) -> Dict[str, str]:
        return dict(tree_file_path=args.tree, verbose=args.verbose)

    def calculate_distance_between_pairs(
        self,
        tips: List[str],
        tree
    ) -> Tuple[
        List[Tuple[str, str]],
        List[float],
    ]:
        combos = list(itertools.combinations(tips, 2))

        patristic_distances = [
            tree.distance(combo[0], combo[1]) for combo in combos
        ]

        return combos, patristic_distances

    def calculate_patristic_distances(
        self,
        tree: Newick.Tree,
    ) -> Tuple[
        List[float],
        List[Tuple[str, str]],
        Dict[str, float],
    ]:
        tips = self.get_tip_names_from_tree(tree)

        combos, patristic_distances = \
            self.calculate_distance_between_pairs(tips, tree)

        stats = calculate_summary_statistics_from_arr(patristic_distances)

        return patristic_distances, combos, stats
