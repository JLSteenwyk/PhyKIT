from __future__ import annotations

import sys
from bisect import bisect_left

from .base import Tree


def calculate_summary_statistics_from_arr(*args, **kwargs):
    from ...helpers.stats_summary import (
        calculate_summary_statistics_from_arr as _calculate_summary_statistics_from_arr,
    )

    return _calculate_summary_statistics_from_arr(*args, **kwargs)


def print_summary_statistics(*args, **kwargs):
    from ...helpers.stats_summary import print_summary_statistics as _print_summary_statistics

    return _print_summary_statistics(*args, **kwargs)


def _json_dumps(*args, **kwargs):
    import json

    return json.dumps(*args, **kwargs)


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()
_NUMPY_THRESHOLD_STATS_MIN_VALUES = 4096
_NUMPY_THRESHOLD_STATS_MIN_THRESHOLDS = 16


class BipartitionSupportStats(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"], verbose=parsed["verbose"])
        self.json_output = parsed["json_output"]
        self.thresholds = parsed["thresholds"]

    def run(self) -> None:
        tree = self.read_tree_file_unmodified()
        if self.verbose:
            bs_vals, term_names = self.get_bipartition_support_vals(tree)
        else:
            bs_vals = self.get_bipartition_support_values(tree)
            term_names = []
        threshold_stats = self.calculate_threshold_stats(bs_vals, self.thresholds)

        if self.json_output:
            payload = self.build_json_output(bs_vals, term_names, threshold_stats)
            try:
                print(_json_dumps(payload, sort_keys=True))
            except BrokenPipeError:
                pass
            return

        if self.verbose:
            try:
                lines = [
                    f"{bs_vals[i]} {';'.join(term_names[i])}"
                    for i in range(len(bs_vals))
                ]
                if lines:
                    print("\n".join(lines))
            except BrokenPipeError:
                pass
        else:
            stats = calculate_summary_statistics_from_arr(bs_vals)
            print_summary_statistics(stats)
            self._print_threshold_stats(threshold_stats)

    def process_args(self, args) -> dict[str, str]:
        return dict(
            tree_file_path=args.tree,
            verbose=args.verbose,
            json_output=args.json,
            thresholds=self.parse_thresholds(args.thresholds),
        )

    def parse_thresholds(self, thresholds_arg: str) -> list[float]:
        if thresholds_arg is None:
            return []

        thresholds: list[float] = []
        for raw_value in thresholds_arg.split(","):
            value = raw_value.strip()
            if not value:
                continue
            try:
                parsed = float(value)
            except ValueError:
                print("All --thresholds values must be numeric.")
                sys.exit(2)
            thresholds.append(parsed)

        if not thresholds:
            print("Provide at least one numeric cutoff with --thresholds.")
            sys.exit(2)

        return thresholds

    @staticmethod
    def _print_threshold_stats(threshold_stats: list[dict[str, float]]) -> None:
        lines = []
        append = lines.append
        for threshold_info in threshold_stats:
            threshold = threshold_info["threshold"]
            count_below = threshold_info["count_below"]
            fraction_below = threshold_info["fraction_below"]
            append(
                f"below {round(threshold, 4)}: {count_below} "
                f"({round(fraction_below * 100, 4)}%)"
            )
        if lines:
            try:
                print("\n".join(lines))
            except BrokenPipeError:
                pass

    def calculate_threshold_stats(
        self, bs_vals: list[float], thresholds: list[float]
    ) -> list[dict[str, float]]:
        if not thresholds:
            return []

        if len(bs_vals) == 0:
            return [
                dict(threshold=threshold, count_below=0, fraction_below=0.0)
                for threshold in thresholds
            ]

        total = len(bs_vals)
        is_numpy_values = type(bs_vals).__module__.startswith("numpy")
        if len(thresholds) == 1:
            threshold = thresholds[0]
            if is_numpy_values:
                count_below = int(np.count_nonzero(bs_vals < threshold))
            elif total >= _NUMPY_THRESHOLD_STATS_MIN_VALUES:
                values = np.asarray(bs_vals, dtype=float)
                count_below = int(np.count_nonzero(values < threshold))
            else:
                count_below = sum(val < threshold for val in bs_vals)
            return [
                dict(
                    threshold=threshold,
                    count_below=count_below,
                    fraction_below=count_below / total,
                )
            ]

        if is_numpy_values:
            sorted_vals = np.sort(bs_vals)
            counts_below = np.searchsorted(
                sorted_vals,
                thresholds,
                side="left",
            )
            return [
                dict(
                    threshold=threshold,
                    count_below=int(count_below),
                    fraction_below=float(count_below / total),
                )
                for threshold, count_below in zip(thresholds, counts_below)
            ]

        if (
            len(thresholds) >= _NUMPY_THRESHOLD_STATS_MIN_THRESHOLDS
            and total >= _NUMPY_THRESHOLD_STATS_MIN_VALUES
        ):
            sorted_vals = np.sort(np.asarray(bs_vals, dtype=float))
            counts_below = np.searchsorted(
                sorted_vals,
                thresholds,
                side="left",
            )
            return [
                dict(
                    threshold=threshold,
                    count_below=int(count_below),
                    fraction_below=float(count_below / total),
                )
                for threshold, count_below in zip(thresholds, counts_below)
            ]

        sorted_vals = sorted(bs_vals)
        return [
            dict(
                threshold=threshold,
                count_below=(count_below := bisect_left(sorted_vals, threshold)),
                fraction_below=count_below / total,
            )
            for threshold in thresholds
        ]

    def build_json_output(
        self,
        bs_vals: list[float],
        term_names: list[list[str]],
        threshold_stats: list[dict[str, float]],
    ) -> dict:
        payload: dict = dict(
            thresholds=threshold_stats,
            verbose=self.verbose,
        )

        if self.verbose:
            payload["bipartitions"] = [
                {"support": support, "terminals": terminals}
                for support, terminals in zip(bs_vals, term_names)
            ]
        else:
            payload["summary"] = calculate_summary_statistics_from_arr(bs_vals)

        return self._to_builtin(payload)

    def _to_builtin(self, value):
        if value is None or type(value) in (str, int, float, bool):
            return value
        if isinstance(value, dict):
            converted = {}
            changed = False
            for key, sub_value in value.items():
                converted_value = self._to_builtin(sub_value)
                changed = changed or converted_value is not sub_value
                converted[key] = converted_value
            return converted if changed else value
        if isinstance(value, list):
            converted = []
            append = converted.append
            changed = False
            for sub_value in value:
                converted_value = self._to_builtin(sub_value)
                changed = changed or converted_value is not sub_value
                append(converted_value)
            return converted if changed else value
        value_type_module = type(value).__module__
        if value_type_module == "numpy" or value_type_module.startswith("numpy."):
            if hasattr(value, "tolist"):
                return self._to_builtin(value.tolist())
            if hasattr(value, "item"):
                return value.item()
        return value

    def get_bipartition_support_vals(
        self,
        tree: Newick.Tree,
    ) -> tuple[list[float], list[list[str]]]:
        direct_result = self._get_bipartition_support_vals_direct(tree)
        if direct_result is not None:
            return direct_result

        # Cache descendant names before the output-order nonterminal traversal.
        bs_vals = []
        term_names = []

        clade_term_names = {}
        for clade in tree.find_clades(order="postorder"):
            if clade.is_terminal():
                clade_term_names[id(clade)] = [clade.name]
            else:
                names = []
                for child in clade.clades:
                    names.extend(clade_term_names.get(id(child), []))
                clade_term_names[id(clade)] = names

        for nonterminal in tree.get_nonterminals():
            if nonterminal.confidence is not None:
                bs_vals.append(nonterminal.confidence)
                term_names.append(clade_term_names.get(id(nonterminal), []))

        return bs_vals, term_names

    def get_bipartition_support_values(self, tree: Newick.Tree):
        direct_result = self._get_bipartition_support_values_array_direct(tree)
        if direct_result is not None:
            return direct_result

        bs_vals, _ = self.get_bipartition_support_vals(tree)
        return bs_vals

    @staticmethod
    def _get_bipartition_support_values_array_direct(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        def iter_support_values():
            stack = [root]
            append = stack.append
            pop = stack.pop
            while stack:
                clade = pop()
                children = clade.clades
                if children:
                    if clade.confidence is not None:
                        yield clade.confidence
                    child_count = len(children)
                    if child_count == 2:
                        append(children[1])
                        append(children[0])
                    else:
                        for idx in range(child_count - 1, -1, -1):
                            append(children[idx])

        try:
            return np.fromiter(iter_support_values(), dtype=float)
        except AttributeError:
            return None

    @staticmethod
    def _get_bipartition_support_vals_direct(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        preorder_internal = []
        preorder = []
        try:
            stack = [root]
            append_stack = stack.append
            pop = stack.pop
            append_preorder = preorder.append
            append_internal = preorder_internal.append
            while stack:
                clade = pop()
                append_preorder(clade)
                children = clade.clades
                if children:
                    append_internal(clade)
                    child_count = len(children)
                    if child_count == 2:
                        append_stack(children[1])
                        append_stack(children[0])
                    else:
                        for idx in range(child_count - 1, -1, -1):
                            append_stack(children[idx])
        except AttributeError:
            return None

        clade_term_names = {}
        try:
            for clade in reversed(preorder):
                children = clade.clades
                if children:
                    names = []
                    for child in children:
                        names.extend(clade_term_names.get(id(child), []))
                    clade_term_names[id(clade)] = names
                else:
                    clade_term_names[id(clade)] = [clade.name]
        except AttributeError:
            return None

        bs_vals = []
        term_names = []
        for nonterminal in preorder_internal:
            if nonterminal.confidence is not None:
                bs_vals.append(nonterminal.confidence)
                term_names.append(clade_term_names.get(id(nonterminal), []))

        return bs_vals, term_names
