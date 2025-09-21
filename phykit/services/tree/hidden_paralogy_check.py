import sys
from typing import Dict, List, Union
from functools import lru_cache
import multiprocessing as mp
from functools import partial

from Bio import Phylo

from .base import Tree


class HiddenParalogyCheck(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    @staticmethod
    def _process_clade_batch(clade_batch, tree_file_path, master_tree_tips):
        """Process a batch of clades in parallel."""
        batch_results = []

        for clade in clade_batch:
            # Read a fresh copy of the tree for each clade
            tree = Phylo.read(tree_file_path, "newick")
            clade_of_interest = set(clade).intersection(master_tree_tips)

            if len(clade_of_interest) <= 1:
                batch_results.append(["insufficient_taxon_representation"])
                continue

            diff_tips = master_tree_tips - clade_of_interest

            # Root and find common ancestor
            try:
                tree.root_with_outgroup(list(diff_tips))
                subtree = tree.common_ancestor(clade_of_interest)

                # Get terminal names efficiently
                common_ancestor_tips = set(tip.name for tip in subtree.get_terminals())

                diff_tips_between_clade_and_curr_tree = \
                    clade_of_interest.symmetric_difference(common_ancestor_tips)

                batch_results.append([
                    "monophyletic" if not diff_tips_between_clade_and_curr_tree else "not_monophyletic",
                    list(diff_tips_between_clade_and_curr_tree),
                ])
            except (ValueError, AttributeError):
                # Handle edge cases where rooting fails
                batch_results.append(["processing_error"])

        return batch_results

    def run(self) -> None:
        # Read the master tree once to get all tip names
        master_tree = self.read_tree_file()
        master_tree_tips = frozenset(self.get_tip_names_from_tree(master_tree))

        # Read clades
        clades = self.read_clades_file(self.clade)

        # For small datasets, process sequentially
        if len(clades) < 10:
            res_arr = []
            for clade in clades:
                # Read a fresh tree for each clade instead of deep copying
                tree = Phylo.read(self.tree_file_path, "newick")
                clade_of_interest = set(clade).intersection(master_tree_tips)

                if len(clade_of_interest) <= 1:
                    res_arr.append(["insufficient_taxon_representation"])
                    continue

                diff_tips = master_tree_tips - clade_of_interest
                tree.root_with_outgroup(list(diff_tips))

                subtree = tree.common_ancestor(clade_of_interest)
                common_ancestor_tips = set(self.get_tip_names_from_tree(subtree))

                diff_tips_between_clade_and_curr_tree = \
                    clade_of_interest.symmetric_difference(common_ancestor_tips)

                res_arr.append([
                    "monophyletic" if not diff_tips_between_clade_and_curr_tree else "not_monophyletic",
                    list(diff_tips_between_clade_and_curr_tree),
                ])
        else:
            # Use multiprocessing for larger datasets
            num_workers = min(mp.cpu_count(), 8)
            batch_size = max(1, len(clades) // num_workers)

            # Create clade batches
            clade_batches = [clades[i:i + batch_size]
                           for i in range(0, len(clades), batch_size)]

            # Process batches in parallel
            process_func = partial(
                self._process_clade_batch,
                tree_file_path=self.tree_file_path,
                master_tree_tips=master_tree_tips
            )

            with mp.Pool(processes=num_workers) as pool:
                batch_results = pool.map(process_func, clade_batches)

            # Flatten results
            res_arr = []
            for batch_result in batch_results:
                res_arr.extend(batch_result)

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
            sys.exit(2)

    def print_results(self, res_arr: List[List[Union[List, str]]]) -> None:
        for res in res_arr:
            print(res[0])
