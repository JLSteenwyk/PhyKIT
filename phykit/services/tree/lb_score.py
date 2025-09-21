import sys
import itertools
from typing import Dict, List, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp
import numpy as np
import pickle
from functools import lru_cache
try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False

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

    @staticmethod
    def _calculate_distances_batch(tree_pickle, tip_pairs):
        """Calculate distances for a batch of tip pairs."""
        tree = pickle.loads(tree_pickle)
        return [tree.distance(tip1, tip2) for tip1, tip2 in tip_pairs]

    def calculate_average_distance_between_tips(
        self,
        tips: List[str],
        tree: Newick.Tree,
    ) -> float:
        num_tips = len(tips)
        if num_tips < 2:
            return 0

        # Get all combinations
        all_pairs = list(itertools.combinations(tips, 2))
        num_combos = len(all_pairs)

        # For small datasets, use sequential processing
        if num_combos < 100:
            total_dist = sum(
                tree.distance(tip1, tip2)
                for tip1, tip2 in all_pairs
            )
        else:
            # Use multiprocessing for large datasets
            tree_pickle = pickle.dumps(tree)
            batch_size = max(50, num_combos // mp.cpu_count())

            with ProcessPoolExecutor(max_workers=min(mp.cpu_count(), 8)) as executor:
                futures = []
                for i in range(0, num_combos, batch_size):
                    batch = all_pairs[i:i + batch_size]
                    futures.append(
                        executor.submit(self._calculate_distances_batch, tree_pickle, batch)
                    )

                total_dist = 0
                # Add progress bar if available and dataset is large
                if HAS_TQDM and num_combos > 1000:
                    futures_iter = tqdm(as_completed(futures), total=len(futures), desc="Computing distances")
                else:
                    futures_iter = as_completed(futures)

                for future in futures_iter:
                    total_dist += sum(future.result())

        return total_dist / num_combos if num_combos else 0

    @staticmethod
    def _calculate_tip_distances_batch(tree_pickle, tips_data):
        """Calculate average distances for a batch of tips."""
        tree = pickle.loads(tree_pickle)
        results = []

        for tip, other_tips in tips_data:
            distances = [tree.distance(tip, other_tip) for other_tip in other_tips]
            avg_dist = sum(distances) / len(distances) if distances else 0
            results.append(avg_dist)

        return results

    def calculate_average_distance_of_taxon_to_other_taxa(
        self,
        tips: List[str],
        tree: Newick.Tree,
    ) -> List[float]:
        # IMPORTANT: Original code has a bug where it uses set(tip) which creates
        # a set of characters, not a set containing the tip. This includes the
        # current tip in distance calculations. We preserve this for compatibility.

        # For small datasets or to maintain exact compatibility, use sequential
        if len(tips) <= 50:
            avg_PDis = []
            for tip in tips:
                # Preserve the original bug: set(tip) creates set of characters
                tips_minus_i = list(set(tips) - set(tip))
                PDi = []
                for tip_minus in tips_minus_i:
                    PDi.append(tree.distance(tip, tip_minus))
                PDi = sum(PDi) / len(PDi) if PDi else 0
                avg_PDis.append(PDi)

            return avg_PDis

        # For larger datasets, use parallel processing but preserve the bug
        tips_data = []
        for tip in tips:
            # Preserve the bug: set(tip) creates set of characters
            tips_minus_i = list(set(tips) - set(tip))
            tips_data.append((tip, tips_minus_i))

        # Process in batches
        batch_size = max(10, len(tips) // mp.cpu_count())
        tree_pickle = pickle.dumps(tree)

        with ProcessPoolExecutor(max_workers=min(mp.cpu_count(), 8)) as executor:
            # Keep track of batch order
            future_to_index = {}

            for i in range(0, len(tips_data), batch_size):
                batch = tips_data[i:i + batch_size]
                future = executor.submit(self._calculate_tip_distances_batch, tree_pickle, batch)
                future_to_index[future] = i

            # Collect results in order
            results_dict = {}
            for future in as_completed(future_to_index):
                batch_index = future_to_index[future]
                results_dict[batch_index] = future.result()

        # Reconstruct ordered results
        avg_PDis = []
        for i in sorted(results_dict.keys()):
            avg_PDis.extend(results_dict[i])

        return avg_PDis

    def calculate_lb_score_per_taxa(
        self,
        avg_PDis: List[float],
        avg_dist: float
    ) -> List[float]:
        if avg_dist == 0:
            try:
                print("Invalid tree. Tree should contain branch lengths")
                sys.exit(2)
            except BrokenPipeError:
                pass
            return []

        # Use NumPy for vectorized computation
        PDis_array = np.array(avg_PDis)
        LBis = ((PDis_array / avg_dist) - 1) * 100

        return LBis.tolist()

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
