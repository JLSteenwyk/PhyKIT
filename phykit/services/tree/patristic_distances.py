from typing import Dict, List, Tuple
import itertools
import multiprocessing as mp
from functools import partial
import pickle
import sys

from Bio.Phylo import Newick
try:
    from tqdm import tqdm
except ImportError:
    # Fallback if tqdm is not installed
    def tqdm(iterable, *args, **kwargs):
        return iterable

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

    def _calculate_distance_batch(self, tree_pickle, combo_batch):
        """Helper function to calculate distances for a batch of combinations."""
        tree = pickle.loads(tree_pickle)
        return [tree.distance(combo[0], combo[1]) for combo in combo_batch]

    def calculate_distance_between_pairs(
        self,
        tips: List[str],
        tree
    ) -> Tuple[
        List[Tuple[str, str]],
        List[float],
    ]:
        combos = list(itertools.combinations(tips, 2))

        # For small datasets, use the original single-threaded approach
        if len(combos) < 100:
            patristic_distances = [
                tree.distance(combo[0], combo[1]) for combo in combos
            ]
        else:
            # Use multiprocessing for larger datasets
            # Serialize the tree once to avoid repeated serialization
            tree_pickle = pickle.dumps(tree)

            # Determine optimal number of workers
            num_workers = min(mp.cpu_count(), 8)

            # Split combos into chunks for parallel processing
            chunk_size = max(1, len(combos) // (num_workers * 4))
            combo_chunks = [combos[i:i + chunk_size] for i in range(0, len(combos), chunk_size)]

            # Create partial function with the pickled tree
            calc_func = partial(self._calculate_distance_batch, tree_pickle)

            # Process in parallel with progress bar
            with mp.Pool(processes=num_workers) as pool:
                # Only show progress bar if stderr is a tty (not redirected)
                if sys.stderr.isatty():
                    results = list(tqdm(
                        pool.imap(calc_func, combo_chunks),
                        total=len(combo_chunks),
                        desc="Calculating patristic distances",
                        unit="batch"
                    ))
                else:
                    results = pool.map(calc_func, combo_chunks)

            # Flatten the results
            patristic_distances = [dist for chunk_result in results for dist in chunk_result]

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
