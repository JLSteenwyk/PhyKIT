from typing import Dict, List, Set, Tuple
from functools import lru_cache
from concurrent.futures import ProcessPoolExecutor
import pickle

from Bio.Phylo import Newick

from .base import Tree


class RobinsonFouldsDistance(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree_zero = self.read_tree_file()
        tree_one = self.read_tree1_file()

        # get shared tree tip names - use sets for efficiency
        tree_zero_tips = set(self.get_tip_names_from_tree(tree_zero))
        tree_one_tips = set(self.get_tip_names_from_tree(tree_one))
        shared_tree_tips = tree_zero_tips & tree_one_tips

        # prune to common set - already have sets
        tree_zero_tips_to_prune = list(tree_zero_tips - shared_tree_tips)
        tree_one_tips_to_prune = list(tree_one_tips - shared_tree_tips)

        if tree_zero_tips_to_prune:
            tree_zero = self.prune_tree_using_taxa_list(tree_zero, tree_zero_tips_to_prune)
        if tree_one_tips_to_prune:
            tree_one = self.prune_tree_using_taxa_list(tree_one, tree_one_tips_to_prune)

        # Get first terminal for rooting
        tip_for_rooting = tree_zero.get_terminals()[0].name
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
        plain_rf = self.compare_trees_optimized(plain_rf, tree_zero, tree_one)
        plain_rf = self.compare_trees_optimized(plain_rf, tree_one, tree_zero)

        tip_count = tree_zero.count_terminals()
        normalized_rf = plain_rf / (2 * (tip_count - 3))

        return plain_rf, normalized_rf

    def compare_trees_optimized(
        self,
        plain_rf: int,
        tree_zero: Newick.Tree,
        tree_one: Newick.Tree
    ) -> int:
        # Cache tip names for clades to avoid recomputation
        tip_names_cache = {}

        def get_cached_tips(clade):
            clade_id = id(clade)
            if clade_id not in tip_names_cache:
                tip_names_cache[clade_id] = frozenset(self.get_tip_names_from_tree(clade))
            return tip_names_cache[clade_id]

        # loop through tree_zero and find similar clade in tree_one
        for clade_zero in tree_zero.get_nonterminals()[1:]:
            # Get tip names from tree_zero clade
            tip_names_zero = get_cached_tips(clade_zero)
            # get common ancestor of tree_zero tip names in tree_one
            clade_one = tree_one.common_ancestor(list(tip_names_zero))
            # Get tip names from tree_one clade
            tip_names_one = get_cached_tips(clade_one)
            # compare the list of tip names
            if tip_names_zero != tip_names_one:
                plain_rf += 1

        return plain_rf

    @staticmethod
    def _calculate_rf_batch(tree_pairs_pickle):
        """Calculate RF distance for a batch of tree pairs in parallel."""
        tree_pairs = pickle.loads(tree_pairs_pickle)
        results = []

        for tree_zero, tree_one in tree_pairs:
            rf_calc = RobinsonFouldsDistance.__new__(RobinsonFouldsDistance)
            rf_calc.__dict__.update({'tree_format': 'newick'})

            # Calculate bipartitions
            bipartitions_zero = rf_calc.get_all_bipartitions(tree_zero)
            bipartitions_one = rf_calc.get_all_bipartitions(tree_one)

            # Calculate RF distance
            plain_rf = len(bipartitions_zero ^ bipartitions_one)  # Symmetric difference
            tip_count = tree_zero.count_terminals()
            normalized_rf = plain_rf / (2 * (tip_count - 3))

            results.append((plain_rf, normalized_rf))

        return results

    def calculate_multiple_rf_distances(self, tree_pairs: List[Tuple]) -> List[Tuple[int, float]]:
        """Calculate RF distances for multiple tree pairs in parallel."""
        if len(tree_pairs) < 5:
            # Sequential for small datasets
            results = []
            for tree_zero, tree_one in tree_pairs:
                plain_rf, normalized_rf = self.calculate_robinson_foulds_distance(tree_zero, tree_one)
                results.append((plain_rf, normalized_rf))
            return results

        # Parallel processing for larger datasets
        batch_size = max(2, len(tree_pairs) // 4)
        batches = [tree_pairs[i:i + batch_size] for i in range(0, len(tree_pairs), batch_size)]

        with ProcessPoolExecutor(max_workers=min(4, len(batches))) as executor:
            futures = []
            for batch in batches:
                batch_pickle = pickle.dumps(batch)
                futures.append(executor.submit(self._calculate_rf_batch, batch_pickle))

            all_results = []
            for future in futures:
                all_results.extend(future.result())

        return all_results
