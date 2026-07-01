from __future__ import annotations

import itertools
import sys

from .base import Tree


def calculate_summary_statistics_from_arr(*args, **kwargs):
    from ...helpers.stats_summary import (
        calculate_summary_statistics_from_arr as _calculate_summary_statistics_from_arr,
    )

    return _calculate_summary_statistics_from_arr(*args, **kwargs)


def print_summary_statistics(*args, **kwargs):
    from ...helpers.stats_summary import print_summary_statistics as _print_summary_statistics

    return _print_summary_statistics(*args, **kwargs)


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class _LazyMultiprocessing:
    def cpu_count(self):
        import multiprocessing as _mp

        return _mp.cpu_count()

    def Pool(self, *args, **kwargs):
        import multiprocessing as _mp

        return _mp.Pool(*args, **kwargs)


mp = _LazyMultiprocessing()


class _LazyPickle:
    def dumps(self, *args, **kwargs):
        import pickle as _pickle

        return _pickle.dumps(*args, **kwargs)

    def loads(self, *args, **kwargs):
        import pickle as _pickle

        return _pickle.loads(*args, **kwargs)


pickle = _LazyPickle()


def _with_optional_progress(iterable, **kwargs):
    try:
        from tqdm import tqdm
    except ImportError:
        return iterable
    return tqdm(iterable, **kwargs)


class PatristicDistances(Tree):
    _PAIRWISE_LCA_DEPTH_THRESHOLD = 64

    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"], verbose=parsed["verbose"])
        self.json_output = parsed["json_output"]

    def run(self):
        tree = self.read_tree_file_unmodified()
        if self.verbose:
            patristic_distances, combos, stats = \
                self.calculate_patristic_distances(tree)
        else:
            patristic_distances = []
            combos = []
            stats = self.calculate_patristic_distance_stats(tree)

        if self.json_output:
            if self.verbose:
                rows = [
                    {
                        "taxon_a": combo[0],
                        "taxon_b": combo[1],
                        "patristic_distance": round(patristic_distance, 4),
                    }
                    for combo, patristic_distance in zip(combos, patristic_distances)
                ]
                print_json(
                    dict(
                        verbose=True,
                        rows=rows,
                        pairs=rows,
                    )
                )
            else:
                print_json(dict(verbose=False, summary=stats))
            return

        if self.verbose:
            try:
                lines = [
                    f"{combo[0]}\t{combo[1]}\t{round(patristic_distance, 4)}"
                    for combo, patristic_distance in zip(
                        combos, patristic_distances
                    )
                ]
                if lines:
                    print("\n".join(lines))
            except BrokenPipeError:
                pass
        else:
            print_summary_statistics(stats)

    def process_args(self, args) -> dict[str, str]:
        return dict(
            tree_file_path=args.tree,
            verbose=args.verbose,
            json_output=getattr(args, "json", False),
        )

    def _calculate_distance_batch(self, tree_pickle, combo_batch):
        """Helper function to calculate distances for a batch of combinations."""
        tree = pickle.loads(tree_pickle)
        return [tree.distance(combo[0], combo[1]) for combo in combo_batch]

    @staticmethod
    def _batched_tip_pairs(tips: list[str], batch_size: int):
        pairs = itertools.combinations(tips, 2)
        while True:
            batch = list(itertools.islice(pairs, batch_size))
            if not batch:
                break
            yield batch

    def calculate_distance_between_pairs(
        self,
        tips: list[str],
        tree
    ) -> tuple[
        list[tuple[str, str]],
        list[float],
    ]:
        fast_result = self.calculate_pairwise_tip_distances_fast(tree, tips)
        if fast_result is not None:
            return fast_result

        combos = list(itertools.combinations(tips, 2))

        # For small datasets, use the original single-threaded approach
        if len(combos) < 100:
            patristic_distances = [
                tree.distance(combo[0], combo[1]) for combo in combos
            ]
        else:
            # Use multiprocessing for larger datasets
            from functools import partial

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
                    results = list(_with_optional_progress(
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

    def calculate_distance_values_between_pairs(
        self,
        tips: list[str],
        tree,
    ) -> list[float]:
        fast_result = self.calculate_pairwise_tip_distance_values_fast(tree, tips)
        if fast_result is not None:
            return fast_result

        return self._calculate_distance_values_between_pairs_fallback(tips, tree)

    def _calculate_distance_values_between_pairs_fallback(
        self,
        tips: list[str],
        tree,
    ) -> list[float]:
        num_pairs = len(tips) * (len(tips) - 1) // 2
        if num_pairs < 100:
            return [
                tree.distance(tip_a, tip_b)
                for tip_a, tip_b in itertools.combinations(tips, 2)
            ]

        from functools import partial

        tree_pickle = pickle.dumps(tree)
        num_workers = min(mp.cpu_count(), 8)
        chunk_size = max(1, num_pairs // (num_workers * 4))
        pair_chunks = self._batched_tip_pairs(tips, chunk_size)
        total_chunks = (num_pairs + chunk_size - 1) // chunk_size
        calc_func = partial(self._calculate_distance_batch, tree_pickle)

        with mp.Pool(processes=num_workers) as pool:
            if sys.stderr.isatty():
                results = list(_with_optional_progress(
                    pool.imap(calc_func, pair_chunks),
                    total=total_chunks,
                    desc="Calculating patristic distances",
                    unit="batch",
                ))
            else:
                results = pool.map(calc_func, pair_chunks)

        return [dist for chunk_result in results for dist in chunk_result]

    def calculate_pairwise_tip_distance_values_fast(
        self,
        tree,
        tips: list[str],
    ):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        parent_indices = []
        levels = []
        depths = []
        terminal_by_name = {}
        stack = [(root, -1, 0, 0.0)]
        try:
            pop = stack.pop
            append = stack.append
            while stack:
                clade, parent_idx, level, depth = pop()
                clade_idx = len(parent_indices)
                parent_indices.append(parent_idx)
                levels.append(level)
                depths.append(depth)
                children = clade.clades
                if children:
                    next_level = level + 1
                    child_count = len(children)
                    if child_count == 2:
                        child = children[1]
                        append(
                            (
                                child,
                                clade_idx,
                                next_level,
                                depth + (child.branch_length or 0.0),
                            )
                        )
                        child = children[0]
                        append(
                            (
                                child,
                                clade_idx,
                                next_level,
                                depth + (child.branch_length or 0.0),
                            )
                        )
                    else:
                        for idx in range(child_count - 1, -1, -1):
                            child = children[idx]
                            append(
                                (
                                    child,
                                    clade_idx,
                                    next_level,
                                    depth + (child.branch_length or 0.0),
                                )
                            )
                else:
                    terminal_by_name[clade.name] = clade_idx
        except (AttributeError, TypeError):
            return None

        try:
            tip_indices = [terminal_by_name[tip] for tip in tips]
        except (KeyError, TypeError):
            return None

        if not tip_indices:
            return []

        max_tip_level = max(levels[tip_idx] for tip_idx in tip_indices)
        if max_tip_level > self._PAIRWISE_LCA_DEPTH_THRESHOLD:
            return self._pairwise_tip_distance_values_from_lca_index(
                tip_indices,
                parent_indices,
                levels,
                depths,
                max_tip_level,
            )

        return self._pairwise_tip_distance_values_from_paths(
            tip_indices,
            parent_indices,
            depths,
        )

    @staticmethod
    def _pairwise_tip_distance_values_from_paths(
        tip_indices,
        parent_indices,
        depths,
    ) -> list[float]:
        tip_paths = []
        for tip_idx in tip_indices:
            path = [tip_idx]
            current = tip_idx
            while parent_indices[current] != -1:
                current = parent_indices[current]
                path.append(current)
            path.reverse()
            tip_paths.append(tuple(path))

        distances = []
        for i in range(len(tip_indices) - 1):
            path_a = tip_paths[i]
            depth_a = depths[tip_indices[i]]
            for j in range(i + 1, len(tip_indices)):
                mrca = 0
                for clade_a, clade_b in zip(path_a, tip_paths[j]):
                    if clade_a != clade_b:
                        break
                    mrca = clade_a
                distances.append(depth_a + depths[tip_indices[j]] - 2 * depths[mrca])

        return distances

    @staticmethod
    def _pairwise_tip_distance_values_from_lca_index(
        tip_indices,
        parent_indices,
        levels,
        depths,
        max_tip_level,
    ) -> list[float]:
        jump_count = max(1, max_tip_level.bit_length())
        ancestors = [parent_indices]
        for _ in range(1, jump_count):
            previous = ancestors[-1]
            ancestors.append([
                previous[parent_idx] if parent_idx != -1 else -1
                for parent_idx in previous
            ])

        def lca_index(node_a, node_b):
            if levels[node_a] < levels[node_b]:
                node_a, node_b = node_b, node_a

            level_diff = levels[node_a] - levels[node_b]
            bit = 0
            while level_diff:
                if level_diff & 1:
                    node_a = ancestors[bit][node_a]
                level_diff >>= 1
                bit += 1

            if node_a == node_b:
                return node_a

            for bit in range(jump_count - 1, -1, -1):
                ancestor_a = ancestors[bit][node_a]
                ancestor_b = ancestors[bit][node_b]
                if ancestor_a != ancestor_b:
                    node_a = ancestor_a
                    node_b = ancestor_b

            return parent_indices[node_a]

        distances = []
        for i in range(len(tip_indices) - 1):
            tip_a = tip_indices[i]
            depth_a = depths[tip_a]
            for j in range(i + 1, len(tip_indices)):
                tip_b = tip_indices[j]
                mrca = lca_index(tip_a, tip_b)
                distances.append(depth_a + depths[tip_b] - 2 * depths[mrca])

        return distances

    def calculate_patristic_distances(
        self,
        tree: Newick.Tree,
    ) -> tuple[
        list[float],
        list[tuple[str, str]],
        dict[str, float],
    ]:
        tips = self.get_tip_names_from_tree(tree)

        combos, patristic_distances = \
            self.calculate_distance_between_pairs(tips, tree)

        stats = calculate_summary_statistics_from_arr(patristic_distances)

        return patristic_distances, combos, stats

    def calculate_patristic_distance_stats(
        self,
        tree: Newick.Tree,
    ) -> dict[str, float]:
        tips = self.get_tip_names_from_tree(tree)
        patristic_distances = self.calculate_distance_values_between_pairs(tips, tree)
        return calculate_summary_statistics_from_arr(patristic_distances)
