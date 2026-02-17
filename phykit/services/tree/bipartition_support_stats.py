import sys
from typing import Dict, List, Tuple
import json

from Bio.Phylo import Newick
import numpy as np

from .base import Tree
from ...helpers.stats_summary import (
    calculate_summary_statistics_from_arr,
    print_summary_statistics,
)


class BipartitionSupportStats(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"], verbose=parsed["verbose"])
        self.json_output = parsed["json_output"]
        self.thresholds = parsed["thresholds"]

    def run(self) -> None:
        tree = self.read_tree_file()
        bs_vals, term_names = self.get_bipartition_support_vals(tree)
        threshold_stats = self.calculate_threshold_stats(bs_vals, self.thresholds)

        if self.json_output:
            payload = self.build_json_output(bs_vals, term_names, threshold_stats)
            try:
                print(json.dumps(payload, sort_keys=True))
            except BrokenPipeError:
                pass
            return

        if self.verbose:
            try:
                for i in range(len(bs_vals)):
                    print(bs_vals[i], ";".join(term_names[i]))
            except BrokenPipeError:
                pass
        else:
            stats = calculate_summary_statistics_from_arr(bs_vals)
            print_summary_statistics(stats)
            for threshold_info in threshold_stats:
                threshold = threshold_info["threshold"]
                count_below = threshold_info["count_below"]
                fraction_below = threshold_info["fraction_below"]
                try:
                    print(
                        f"below {round(threshold, 4)}: {count_below} "
                        f"({round(fraction_below * 100, 4)}%)"
                    )
                except BrokenPipeError:
                    pass

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            tree_file_path=args.tree,
            verbose=args.verbose,
            json_output=args.json,
            thresholds=self.parse_thresholds(args.thresholds),
        )

    def parse_thresholds(self, thresholds_arg: str) -> List[float]:
        if thresholds_arg is None:
            return []

        thresholds: List[float] = []
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

    def calculate_threshold_stats(
        self, bs_vals: List[float], thresholds: List[float]
    ) -> List[Dict[str, float]]:
        if len(bs_vals) == 0:
            return [
                dict(threshold=threshold, count_below=0, fraction_below=0.0)
                for threshold in thresholds
            ]

        total = len(bs_vals)
        return [
            dict(
                threshold=threshold,
                count_below=sum(val < threshold for val in bs_vals),
                fraction_below=sum(val < threshold for val in bs_vals) / total,
            )
            for threshold in thresholds
        ]

    def build_json_output(
        self,
        bs_vals: List[float],
        term_names: List[List[str]],
        threshold_stats: List[Dict[str, float]],
    ) -> Dict:
        payload: Dict = dict(
            thresholds=threshold_stats,
            verbose=self.verbose,
        )

        if self.verbose:
            payload["bipartitions"] = [
                dict(support=bs_vals[i], terminals=term_names[i])
                for i in range(len(bs_vals))
            ]
        else:
            payload["summary"] = calculate_summary_statistics_from_arr(bs_vals)

        return self._to_builtin(payload)

    def _to_builtin(self, value):
        if isinstance(value, dict):
            return {k: self._to_builtin(v) for k, v in value.items()}
        if isinstance(value, list):
            return [self._to_builtin(v) for v in value]
        if isinstance(value, np.integer):
            return int(value)
        if isinstance(value, np.floating):
            return float(value)
        return value

    def get_bipartition_support_vals(
        self,
        tree: Newick.Tree,
    ) -> Tuple[List[float], List[List[str]]]:
        # Single pass through nonterminals to avoid duplicate tree traversal
        bs_vals = []
        term_names = []

        # Cache terminals for each nonterminal in one pass
        for nonterminal in tree.get_nonterminals():
            if nonterminal.confidence is not None:
                bs_vals.append(nonterminal.confidence)
                # Get terminal names once for this nonterminal
                term_names.append([term.name for term in nonterminal.get_terminals()])

        return bs_vals, term_names
