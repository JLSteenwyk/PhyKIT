from __future__ import annotations

import math

from .base import Tree


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class DVMC(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.json_output = parsed["json_output"]

    def run(self):
        summary = self._get_simple_newick_terminal_distance_stats(
            self.tree_file_path,
            "tree_file_path",
        )
        if summary is None:
            tree = self.read_tree_file_unmodified()
            dvmc = self.determine_dvmc(tree)
        else:
            count, total, total_sq = summary
            avg_dist = total / count
            squared_diff_sum = total_sq - count * (avg_dist * avg_dist)
            dvmc = math.sqrt(squared_diff_sum / (count - 1))
        if self.json_output:
            print_json(dict(dvmc=round(dvmc, 4)))
            return
        print(round(dvmc, 4))

    def process_args(self, args) -> dict[str, str]:
        return dict(tree_file_path=args.tree, json_output=getattr(args, "json", False))

    def determine_dvmc(self, tree: Newick.Tree) -> float:
        dvmc = self._determine_dvmc_standard_tree(tree)
        if dvmc is not None:
            return dvmc

        distances = self.calculate_terminal_root_distances_fast(tree)
        if distances is None:
            terminals = tree.get_terminals()
            try:
                depths = tree.depths()
                root_depth = depths[tree.root]
                distances = [depths[term] - root_depth for term in terminals]
            except (AttributeError, KeyError, TypeError):
                distances = [tree.distance(term) for term in terminals]

        num_spp = len(distances)
        total = 0.0
        total_sq = 0.0
        for distance in distances:
            total += distance
            total_sq += distance * distance
        avg_dist = total / num_spp

        squared_diff_sum = total_sq - num_spp * (avg_dist * avg_dist)

        return math.sqrt(squared_diff_sum / (num_spp - 1))

    @staticmethod
    def _determine_dvmc_standard_tree(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        count = 0
        total = 0.0
        total_sq = 0.0
        stack = [(root, 0.0)]
        pop = stack.pop
        append = stack.append
        try:
            while stack:
                clade, distance = pop()
                children = clade.clades
                if children:
                    for child in children:
                        append(
                            (
                                child,
                                distance + (child.branch_length or 0.0),
                            )
                        )
                else:
                    count += 1
                    total += distance
                    total_sq += distance * distance
        except AttributeError:
            return None

        if count < 2:
            return None

        avg_dist = total / count
        squared_diff_sum = total_sq - count * (avg_dist * avg_dist)
        return math.sqrt(squared_diff_sum / (count - 1))
