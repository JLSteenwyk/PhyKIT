from __future__ import annotations

import sys
import itertools

from .base import Tree


class _LazyPickle:
    _module = None

    def _load(self):
        module = self._module
        if module is None:
            import pickle as module

            self._module = module
        return module

    def dumps(self, *args, **kwargs):
        module = self._load()
        dumps = module.dumps
        self.dumps = dumps
        if "loads" not in self.__dict__:
            self.loads = module.loads

        return dumps(*args, **kwargs)

    def loads(self, *args, **kwargs):
        module = self._load()
        loads = module.loads
        self.loads = loads
        if "dumps" not in self.__dict__:
            self.dumps = module.dumps

        return loads(*args, **kwargs)


pickle = _LazyPickle()
_DENOMINATOR_SCAN_MIN_TIPS = 10_000


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
    _module = None

    def _load(self):
        module = self._module
        if module is None:
            import multiprocessing as module

            self._module = module
        return module

    def cpu_count(self):
        return self._load().cpu_count()


mp = _LazyMultiprocessing()


def ProcessPoolExecutor(*args, **kwargs):
    from concurrent.futures import ProcessPoolExecutor as _ProcessPoolExecutor

    return _ProcessPoolExecutor(*args, **kwargs)


def as_completed(*args, **kwargs):
    from concurrent.futures import as_completed as _as_completed

    return _as_completed(*args, **kwargs)


def _as_completed_with_optional_progress(futures, num_combos):
    futures_iter = as_completed(futures)
    if num_combos <= 1000:
        return futures_iter

    try:
        from tqdm import tqdm
    except ImportError:
        return futures_iter

    return tqdm(
        futures_iter,
        total=len(futures),
        desc="Computing distances",
    )


class LBScore(Tree):
    MP_MIN_DISTANCE_PAIRS = 500_000
    MAX_MP_WORKERS = 8

    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"], verbose=parsed["verbose"])
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        tree = self.read_tree_file_unmodified()
        tips, LBis = self.calculate_lb_score(tree)
        if self.json_output:
            if self.verbose:
                rows = [
                    {"taxon": tip, "lb_score": round(LBi, 4)}
                    for tip, LBi in zip(tips, LBis)
                ]
                print_json(
                    dict(
                        verbose=True,
                        rows=rows,
                        taxa=rows,
                    )
                )
            else:
                stats = calculate_summary_statistics_from_arr(LBis)
                print_json(dict(verbose=False, summary=stats))
            return

        if self.verbose:
            try:
                lines = [
                    f"{tip}\t{round(LBi, 4)}"
                    for tip, LBi in zip(tips, LBis)
                ]
                if lines:
                    print("\n".join(lines))
            except BrokenPipeError:
                pass
        else:
            stats = calculate_summary_statistics_from_arr(LBis)
            print_summary_statistics(stats)

    def process_args(self, args) -> dict[str, str]:
        return dict(
            tree_file_path=args.tree,
            verbose=args.verbose,
            json_output=getattr(args, "json", False),
        )

    @staticmethod
    def _calculate_distances_batch(tree_pickle, tip_pairs):
        """Calculate distances for a batch of tip pairs."""
        tree = pickle.loads(tree_pickle)
        return [tree.distance(tip1, tip2) for tip1, tip2 in tip_pairs]

    @staticmethod
    def _batched_tip_pairs(tips: list[str], batch_size: int):
        pairs = itertools.combinations(tips, 2)
        while True:
            batch = list(itertools.islice(pairs, batch_size))
            if not batch:
                break
            yield batch

    def calculate_average_distance_between_tips(
        self,
        tips: list[str],
        tree: Newick.Tree,
        pairwise_distances: (
            tuple[list[tuple[str, str]], list[float]] | None
        ) = None,
    ) -> float:
        num_tips = len(tips)
        if num_tips < 2:
            return 0

        num_combos = num_tips * (num_tips - 1) // 2
        fast_result = (
            pairwise_distances
            if pairwise_distances is not None
            else self.calculate_pairwise_tip_distances_fast(tree, tips)
        )
        if fast_result is not None:
            _, distances = fast_result
            return sum(distances) / num_combos if num_combos else 0

        # For small/medium fallback datasets, multiprocessing overhead dominates.
        if num_combos < self.MP_MIN_DISTANCE_PAIRS:
            total_dist = sum(
                tree.distance(tip1, tip2)
                for tip1, tip2 in itertools.combinations(tips, 2)
            )
        else:
            # Use multiprocessing for large datasets
            tree_pickle = pickle.dumps(tree)
            cpu_count = mp.cpu_count()
            batch_size = max(50, num_combos // cpu_count)

            with ProcessPoolExecutor(max_workers=min(cpu_count, self.MAX_MP_WORKERS)) as executor:
                futures = []
                for batch in self._batched_tip_pairs(tips, batch_size):
                    futures.append(
                        executor.submit(self._calculate_distances_batch, tree_pickle, batch)
                    )

                total_dist = 0
                for future in _as_completed_with_optional_progress(futures, num_combos):
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
        tips: list[str],
        tree: Newick.Tree,
        pairwise_distances: (
            tuple[list[tuple[str, str]], list[float]] | None
        ) = None,
    ) -> list[float]:
        # IMPORTANT: Original code has a bug where it uses set(tip) which creates
        # a set of characters, not a set containing the tip. This includes the
        # current tip in distance calculations. We preserve this for compatibility.

        fast_result = (
            pairwise_distances
            if pairwise_distances is not None
            else self.calculate_pairwise_tip_distances_fast(tree, tips)
        )
        if fast_result is not None:
            combos, distances = fast_result
            distance_sums = {tip: 0.0 for tip in tips}
            for (tip_a, tip_b), distance in zip(combos, distances):
                distance_sums[tip_a] += distance
                distance_sums[tip_b] += distance

            tip_set = set(tips)
            tip_count = len(tips)
            avg_PDis = []
            for tip in tips:
                # Preserve the original bug: set(tip) creates set of characters
                denominator = self._historical_other_taxa_denominator(
                    tip,
                    tip_set,
                    tip_count,
                )
                avg_PDis.append(
                    distance_sums[tip] / denominator if denominator else 0
                )
            return avg_PDis

        # For small datasets or to maintain exact compatibility, use sequential
        tip_set = set(tips)
        if len(tips) <= 50:
            avg_PDis = []
            for tip in tips:
                # Preserve the original bug: set(tip) creates set of characters
                tips_minus_i = list(tip_set - set(tip))
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
            tips_minus_i = list(tip_set - set(tip))
            tips_data.append((tip, tips_minus_i))

        # Process in batches
        cpu_count = mp.cpu_count()
        batch_size = max(10, len(tips) // cpu_count)
        tree_pickle = pickle.dumps(tree)

        with ProcessPoolExecutor(max_workers=min(cpu_count, self.MAX_MP_WORKERS)) as executor:
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
        avg_PDis: list[float],
        avg_dist: float
    ) -> list[float]:
        if avg_dist == 0:
            try:
                print("Invalid tree. Tree should contain branch lengths")
                sys.exit(2)
            except BrokenPipeError:
                pass
            return []

        scale = 100.0 / avg_dist
        return [(avg_PDi * scale) - 100.0 for avg_PDi in avg_PDis]

    @staticmethod
    def _historical_other_taxa_denominator(
        tip: str,
        tip_set: set[str],
        tip_count: int,
    ) -> int:
        if tip_count < _DENOMINATOR_SCAN_MIN_TIPS:
            return len(tip_set - set(tip))

        denominator = tip_count
        matched = None
        for character in tip:
            if character in tip_set and (
                matched is None or character not in matched
            ):
                denominator -= 1
                if matched is None:
                    matched = {character}
                else:
                    matched.add(character)
        return denominator

    @staticmethod
    def _calculate_lb_components_fast(
        tree: Newick.Tree,
        tips: list[str],
    ) -> tuple[float, list[float]] | None:
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        tip_count = len(tips)
        if tip_count < 2:
            return None

        tip_set = set(tips)
        if len(tip_set) != tip_count:
            return None
        subtree_tip_counts = {}
        subtree_distance_sums = {}
        branch_lengths = {root: 0.0}
        terminal_names_by_node = {}
        stack = [(root, False)]

        try:
            pop = stack.pop
            append = stack.append
            while stack:
                node, visited = pop()
                children = node.clades
                if not visited:
                    append((node, True))
                    child_count = len(children)
                    if child_count == 2:
                        child = children[1]
                        branch_lengths[child] = child.branch_length or 0.0
                        append((child, False))
                        child = children[0]
                        branch_lengths[child] = child.branch_length or 0.0
                        append((child, False))
                    elif child_count:
                        for idx in range(child_count - 1, -1, -1):
                            child = children[idx]
                            branch_lengths[child] = child.branch_length or 0.0
                            append((child, False))
                    continue

                if children:
                    child_tip_count = 0
                    child_distance_sum = 0.0
                    for child in children:
                        child_tip_count += subtree_tip_counts[child]
                        child_distance_sum += (
                            subtree_distance_sums[child]
                            + branch_lengths[child] * subtree_tip_counts[child]
                        )
                    subtree_tip_counts[node] = child_tip_count
                    subtree_distance_sums[node] = child_distance_sum
                elif node.name in tip_set:
                    subtree_tip_counts[node] = 1
                    subtree_distance_sums[node] = 0.0
                    terminal_names_by_node[node] = node.name
                else:
                    subtree_tip_counts[node] = 0
                    subtree_distance_sums[node] = 0.0
        except (AttributeError, TypeError):
            return None

        if subtree_tip_counts[root] != tip_count:
            return None

        total_pairwise_distance = 0.0
        distance_sums_by_node = {root: subtree_distance_sums[root]}
        stack = [root]
        try:
            while stack:
                node = stack.pop()
                for child in node.clades:
                    child_tip_count = subtree_tip_counts[child]
                    branch_length = branch_lengths[child]
                    total_pairwise_distance += (
                        branch_length
                        * child_tip_count
                        * (tip_count - child_tip_count)
                    )
                    distance_sums_by_node[child] = (
                        distance_sums_by_node[node]
                        + branch_length * (tip_count - 2 * child_tip_count)
                    )
                    stack.append(child)
        except (AttributeError, TypeError):
            return None

        avg_dist = total_pairwise_distance / (tip_count * (tip_count - 1) / 2)
        distance_sums_by_name = {
            name: distance_sums_by_node[node]
            for node, name in terminal_names_by_node.items()
        }
        avg_PDis = []
        for tip in tips:
            # Preserve the current optimized path's compatibility behavior:
            # the numerator includes all pairwise distances for the tip, while
            # the denominator retains the historical set(tip) quirk.
            denominator = LBScore._historical_other_taxa_denominator(
                tip,
                tip_set,
                tip_count,
            )
            avg_PDis.append(
                distance_sums_by_name[tip] / denominator if denominator else 0
            )

        return avg_dist, avg_PDis

    def calculate_lb_score(
        self,
        tree: Newick.Tree
    ) -> tuple[list[str], list[float]]:
        tips = self.get_tip_names_from_tree(tree)
        fast_components = self._calculate_lb_components_fast(tree, tips)
        if fast_components is not None:
            avg_dist, avg_PDis = fast_components
            LBis = self.calculate_lb_score_per_taxa(avg_PDis, avg_dist)
            return tips, LBis

        pairwise_distances = self.calculate_pairwise_tip_distances_fast(tree, tips)

        avg_dist = self.calculate_average_distance_between_tips(
            tips,
            tree,
            pairwise_distances,
        )

        avg_PDis = self.calculate_average_distance_of_taxon_to_other_taxa(
            tips,
            tree,
            pairwise_distances,
        )

        LBis = self.calculate_lb_score_per_taxa(avg_PDis, avg_dist)

        return tips, LBis
