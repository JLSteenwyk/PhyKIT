import copy
import numpy as np
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from functools import lru_cache
import pickle

from scipy.stats import pearsonr, zscore

from .base import Tree


class CovaryingEvolutionaryRates(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree_zero = self.read_tree_file()
        tree_one = self.read_tree1_file()
        tree_ref = self.read_reference_tree_file()

        # - Calculate correlation between two gene trees
        # and save results to an array, corrArr.
        # - Branch lengths will also be part of output

        # get tree tip names
        tree_zero_tips = self.get_tip_names_from_tree(tree_zero)
        tree_one_tips = self.get_tip_names_from_tree(tree_one)
        tree_ref_tips = self.get_tip_names_from_tree(tree_ref)

        # get shared tips between the two trees
        shared_tree_tips = self.shared_tips(tree_zero_tips, tree_one_tips)

        # find differences between tree tips and shared tips
        # to determine what tips to prune
        tree_zero_tips_to_prune = list(set(tree_zero_tips) - set(shared_tree_tips))
        tree_one_tips_to_prune = list(set(tree_one_tips) - set(shared_tree_tips))
        tree_ref_tips_to_prune = list(set(tree_ref_tips) - set(shared_tree_tips))

        # get a set of pruned trees
        tree_zero = self.prune_tips(tree_zero, tree_zero_tips_to_prune)
        tree_one = self.prune_tips(tree_one, tree_one_tips_to_prune)
        tree_ref = self.prune_tips(tree_ref, tree_ref_tips_to_prune)

        # obtain corrected branch lengths where branch lengths
        # are corrected by the species tree branch length
        (
            tree_zero_corr_branch_lengths,
            tree_one_corr_branch_lengths,
            tip_names,
        ) = self.correct_branch_lengths(tree_zero, tree_one, tree_ref)

        # remove corrected BLs greater than 5
        outlier_indices = []
        outlier_indices = self.get_indices_of_outlier_branch_lengths(
            tree_zero_corr_branch_lengths, outlier_indices
        )
        outlier_indices = self.get_indices_of_outlier_branch_lengths(
            tree_one_corr_branch_lengths, outlier_indices
        )

        tree_zero_corr_branch_lengths = self.remove_outliers_based_on_indices(
            tree_zero_corr_branch_lengths, outlier_indices
        )
        tree_one_corr_branch_lengths = self.remove_outliers_based_on_indices(
            tree_one_corr_branch_lengths, outlier_indices
        )
        tip_names = self.remove_outliers_based_on_indices(tip_names, outlier_indices)

        # standardize values for final correction
        tree_zero_corr_branch_lengths = zscore(tree_zero_corr_branch_lengths)
        tree_one_corr_branch_lengths = zscore(tree_one_corr_branch_lengths)

        # Calculate correlation and append to results array
        # also keep a list of p values
        corr = list(
            pearsonr(tree_zero_corr_branch_lengths, tree_one_corr_branch_lengths)
        )

        try:
            if self.verbose:
                for val_zero, val_one, tip_name in zip(
                    tree_zero_corr_branch_lengths,
                    tree_one_corr_branch_lengths,
                    tip_names,
                ):
                    print(
                        f"{round(val_zero, 4)}\t{round(val_one, 4)}\t{';'.join(tip_name)}"
                    )
            else:
                print(f"{round(corr[0], 4)}\t{round(corr[1], 6)}")
        except BrokenPipeError:
            pass

    def process_args(self, args):
        return dict(
            tree_file_path=args.tree_zero,
            tree1_file_path=args.tree_one,
            reference=args.reference,
            verbose=args.verbose,
        )

    def get_indices_of_outlier_branch_lengths(
        self, corr_branch_lengths, outlier_indices
    ):
        """
        create index for branch lengths that
        have an absolute value greater than 5
        """
        # Convert to numpy array for vectorized operations
        arr = np.array(corr_branch_lengths, dtype=float)

        # Find outliers using vectorized operations
        new_outliers = np.where((np.abs(arr) > 5) | np.isnan(arr))[0]

        # Combine with existing outliers
        all_outliers = set(outlier_indices)
        all_outliers.update(new_outliers.tolist())

        return list(all_outliers)

    def remove_outliers_based_on_indices(self, corr_branch_lengths, outlier_indices):
        """
        remove value if the value is an outlier according
        to the outlier indices list
        """
        if not outlier_indices:
            return corr_branch_lengths

        # Use numpy for efficient filtering
        mask = np.ones(len(corr_branch_lengths), dtype=bool)
        mask[list(outlier_indices)] = False

        if isinstance(corr_branch_lengths[0], (list, tuple)):
            # Handle list of lists (tip_names)
            return [item for i, item in enumerate(corr_branch_lengths) if mask[i]]
        else:
            # Handle numeric lists
            return [item for i, item in enumerate(corr_branch_lengths) if mask[i]]

    def prune_tips(self, tree, tips):
        """
        prune tips from trees
        """

        for tip in tips:
            tree.prune(tip)

        return tree

    @staticmethod
    def _process_terminal_batch(tree0_pickle, tree1_pickle, terminals_data):
        """Process a batch of terminals in parallel."""
        t0 = pickle.loads(tree0_pickle)
        t1 = pickle.loads(tree1_pickle)

        results = []
        for terminal_name, terminal_bl, sp_tips in terminals_data:
            try:
                newtree = t0.common_ancestor(terminal_name)
                newtree1 = t1.common_ancestor(terminal_name)

                bl0 = round(newtree.branch_length / terminal_bl, 6) if newtree.branch_length else None
                bl1 = round(newtree1.branch_length / terminal_bl, 6) if newtree1.branch_length else None

                if bl0 is not None and bl1 is not None:
                    results.append((bl0, bl1, sp_tips))
            except:
                continue
        return results

    @staticmethod
    def _process_nonterminal_batch(tree0_pickle, tree1_pickle, nonterminals_data):
        """Process a batch of nonterminals in parallel."""
        t0 = pickle.loads(tree0_pickle)
        t1 = pickle.loads(tree1_pickle)

        results = []
        for sp_tips, nonterminal_bl in nonterminals_data:
            try:
                newtree = t0.common_ancestor(sp_tips)
                newtree1 = t1.common_ancestor(sp_tips)

                if newtree.branch_length and newtree1.branch_length and nonterminal_bl:
                    bl0 = round(newtree.branch_length / nonterminal_bl, 6)
                    bl1 = round(newtree1.branch_length / nonterminal_bl, 6)
                    results.append((bl0, bl1, sp_tips))
            except:
                continue
        return results

    def correct_branch_lengths(self, t0, t1, sp):
        """
        obtain a list of corrected branch lengths with parallel processing
        """
        l0 = []
        l1 = []
        tip_names = []

        # Collect terminal data
        terminals = sp.get_terminals()
        nonterminals = sp.get_nonterminals()

        # Process sequentially if small dataset or use parallel processing
        if len(terminals) + len(nonterminals) < 50:
            # Original sequential processing for small datasets
            for i in terminals:
                sp_tips = self.get_tip_names_from_tree(i)
                tip_names.append(sp_tips)
                try:
                    newtree = t0.common_ancestor(i.name)
                    newtree1 = t1.common_ancestor(i.name)
                    if newtree.branch_length and i.branch_length:
                        l0.append(round(newtree.branch_length / i.branch_length, 6))
                        l1.append(round(newtree1.branch_length / i.branch_length, 6))
                except:
                    continue

            for i in nonterminals:
                sp_tips = self.get_tip_names_from_tree(i)
                try:
                    newtree = t0.common_ancestor(sp_tips)
                    newtree1 = t1.common_ancestor(sp_tips)
                    if newtree.branch_length and newtree1.branch_length and i.branch_length:
                        l0.append(round(newtree.branch_length / i.branch_length, 6))
                        l1.append(round(newtree1.branch_length / i.branch_length, 6))
                        tip_names.append(sp_tips)
                except:
                    continue
        else:
            # Parallel processing for large datasets
            tree0_pickle = pickle.dumps(t0)
            tree1_pickle = pickle.dumps(t1)

            # Prepare terminal data
            terminals_data = []
            for i in terminals:
                sp_tips = self.get_tip_names_from_tree(i)
                if i.branch_length:
                    terminals_data.append((i.name, i.branch_length, sp_tips))

            # Prepare nonterminal data
            nonterminals_data = []
            for i in nonterminals:
                if i.branch_length:
                    sp_tips = self.get_tip_names_from_tree(i)
                    nonterminals_data.append((sp_tips, i.branch_length))

            # Process in batches
            batch_size = max(10, (len(terminals_data) + len(nonterminals_data)) // 4)

            try:
                with ProcessPoolExecutor(max_workers=min(4, len(terminals_data) + len(nonterminals_data) // 10)) as executor:
                    futures = []

                    # Submit terminal batches
                    for i in range(0, len(terminals_data), batch_size):
                        batch = terminals_data[i:i+batch_size]
                        futures.append(executor.submit(self._process_terminal_batch, tree0_pickle, tree1_pickle, batch))

                    # Submit nonterminal batches
                    for i in range(0, len(nonterminals_data), batch_size):
                        batch = nonterminals_data[i:i+batch_size]
                        futures.append(executor.submit(self._process_nonterminal_batch, tree0_pickle, tree1_pickle, batch))

                    for future in as_completed(futures):
                        batch_results = future.result()
                        for bl0, bl1, sp_tips in batch_results:
                            l0.append(bl0)
                            l1.append(bl1)
                            tip_names.append(sp_tips)
            except (OSError, ValueError, RuntimeError):
                for i in terminals:
                    sp_tips = self.get_tip_names_from_tree(i)
                    tip_names.append(sp_tips)
                    try:
                        newtree = t0.common_ancestor(i.name)
                        newtree1 = t1.common_ancestor(i.name)
                        if newtree.branch_length and i.branch_length:
                            l0.append(round(newtree.branch_length / i.branch_length, 6))
                            l1.append(round(newtree1.branch_length / i.branch_length, 6))
                    except Exception:
                        continue

                for i in nonterminals:
                    sp_tips = self.get_tip_names_from_tree(i)
                    try:
                        newtree = t0.common_ancestor(sp_tips)
                        newtree1 = t1.common_ancestor(sp_tips)
                        if newtree.branch_length and newtree1.branch_length and i.branch_length:
                            l0.append(round(newtree.branch_length / i.branch_length, 6))
                            l1.append(round(newtree1.branch_length / i.branch_length, 6))
                            tip_names.append(sp_tips)
                    except Exception:
                        continue

        return (l0, l1, tip_names)
