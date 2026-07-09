from __future__ import annotations

import sys

from .base import Tree


def calculate_summary_statistics_from_arr(*args, **kwargs):
    from ...helpers.stats_summary import (
        calculate_summary_statistics_from_arr as _calculate_summary_statistics_from_arr,
    )

    return _calculate_summary_statistics_from_arr(*args, **kwargs)


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def read_single_column_file_to_list(*args, **kwargs):
    from ...helpers.files import (
        read_single_column_file_to_list as _read_single_column_file_to_list,
    )

    return _read_single_column_file_to_list(*args, **kwargs)


_DIRECT_TRAVERSAL_FAILED = object()


class MonophylyCheck(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"], list_of_taxa=parsed["list_of_taxa"])
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        tree = self.read_tree_file_unmodified()
        taxa = read_single_column_file_to_list(self.list_of_taxa)

        res_arr = []

        # Use frozenset for more efficient set operations
        tree_tips = frozenset(self.get_tip_names_from_tree(tree))
        taxa_set = frozenset(taxa)
        taxa_of_interest = taxa_set.intersection(tree_tips)

        if len(taxa_of_interest) <= 1:
            res_arr.append(["insufficient_taxon_representation"])
            sys.exit(2)

        tree, common_ancestor_tips = self._resolve_interest_clade(
            tree, taxa_of_interest, tree_tips
        )
        diff_tips_between_clade_and_curr_tree = list(
            taxa_of_interest.symmetric_difference(common_ancestor_tips)
        )

        stats = self.get_bootstrap_statistics(tree)

        res_arr = self.populate_res_arr(
            diff_tips_between_clade_and_curr_tree, stats, res_arr
        )

        self.print_results(res_arr)

    def process_args(self, args) -> dict[str, str]:
        return dict(
            tree_file_path=args.tree,
            list_of_taxa=args.list_of_taxa,
            json_output=getattr(args, "json", False),
        )

    def _resolve_interest_clade(self, tree, taxa_of_interest, tree_tips):
        if taxa_of_interest == tree_tips:
            try:
                return tree.root, taxa_of_interest
            except AttributeError:
                pass

        clade = self._find_exact_clade_by_taxa(tree, taxa_of_interest)
        if clade is not None:
            return clade, taxa_of_interest

        if taxa_of_interest.issubset(tree_tips):
            shared_tree_tips = list(taxa_of_interest)
            diff_tips = list(tree_tips.difference(taxa_of_interest))
        else:
            shared_tree_tips = self.shared_tips(list(taxa_of_interest), list(tree_tips))
            diff_tips = list(tree_tips - frozenset(shared_tree_tips))

        mutable_tree = self._fast_copy(tree)
        mutable_tree.root_with_outgroup(diff_tips)
        clade = mutable_tree.common_ancestor(shared_tree_tips)
        common_ancestor_tips = frozenset(self.get_tip_names_from_tree(clade))
        return clade, common_ancestor_tips

    @staticmethod
    def _find_exact_clade_by_taxa(tree, taxa):
        target = frozenset(taxa)
        target_size = len(target)
        if target_size == 0:
            return None

        direct_result = MonophylyCheck._find_exact_clade_by_taxa_direct(
            tree, target, target_size
        )
        if direct_result is not _DIRECT_TRAVERSAL_FAILED:
            return direct_result

        selected_counts = {}
        terminal_counts = {}
        seen_terminal_names = set()
        try:
            for clade in tree.find_clades(order="postorder"):
                children = clade.clades
                if children:
                    selected_count = 0
                    terminal_count = 0
                    for child in children:
                        child_id = id(child)
                        selected_count += selected_counts[child_id]
                        terminal_count += terminal_counts[child_id]
                else:
                    name = clade.name
                    if name in seen_terminal_names:
                        return MonophylyCheck._find_exact_clade_by_taxa_sets(
                            tree, target
                        )
                    seen_terminal_names.add(name)
                    selected_count = 1 if name in target else 0
                    terminal_count = 1

                if selected_count == target_size and terminal_count == target_size:
                    return clade

                clade_id = id(clade)
                selected_counts[clade_id] = selected_count
                terminal_counts[clade_id] = terminal_count
        except (AttributeError, KeyError, TypeError):
            return None
        return None

    @staticmethod
    def _find_exact_clade_by_taxa_direct(tree, target, target_size):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return _DIRECT_TRAVERSAL_FAILED

        seen_terminal_names = set()
        stack = [[root, 0, 0, 0]]
        try:
            while stack:
                frame = stack[-1]
                clade = frame[0]
                children = clade.clades

                if children and frame[1] < len(children):
                    child = children[frame[1]]
                    frame[1] += 1
                    stack.append([child, 0, 0, 0])
                    continue

                if children:
                    selected_count = frame[2]
                    terminal_count = frame[3]
                else:
                    name = clade.name
                    if name in seen_terminal_names:
                        return MonophylyCheck._find_exact_clade_by_taxa_sets(
                            tree, target
                        )
                    seen_terminal_names.add(name)
                    selected_count = 1 if name in target else 0
                    terminal_count = 1

                if selected_count == target_size and terminal_count == target_size:
                    return clade

                stack.pop()
                if stack:
                    parent = stack[-1]
                    parent[2] += selected_count
                    parent[3] += terminal_count
        except (AttributeError, KeyError, TypeError):
            return _DIRECT_TRAVERSAL_FAILED
        return None

    @staticmethod
    def _find_exact_clade_by_taxa_sets(tree, taxa):
        try:
            target = frozenset(taxa)
            clade_tips = {}
            for clade in tree.find_clades(order="postorder"):
                if clade.is_terminal():
                    tips = frozenset([clade.name])
                else:
                    tips = frozenset().union(
                        *(clade_tips[id(child)] for child in clade.clades)
                    )
                clade_tips[id(clade)] = tips
                if tips == target:
                    return clade
        except (AttributeError, KeyError, TypeError):
            return None
        return None

    def get_bootstrap_statistics(
        self,
        clade: Newick.Clade
    ) -> dict[str, int | float]:
        bs_vals = self._collect_bootstrap_values_direct(clade)
        if bs_vals is None:
            bs_vals = [
                terminal.confidence for terminal in clade.get_nonterminals()
                if terminal.confidence is not None
            ]

        return calculate_summary_statistics_from_arr(bs_vals)

    @staticmethod
    def _collect_bootstrap_values_direct(clade):
        try:
            clade.clades
        except AttributeError:
            return None

        bs_vals = []
        stack = [clade]
        pop = stack.pop
        extend = stack.extend
        append = bs_vals.append
        try:
            while stack:
                node = pop()
                children = node.clades
                if children:
                    confidence = node.confidence
                    if confidence is not None:
                        append(confidence)
                    extend(children)
        except AttributeError:
            return None

        return bs_vals

    def populate_res_arr(
        self,
        diff_tips_between_clade_and_curr_tree: list[str],
        stats: dict[str, float],
        res_arr: list,
    ) -> list[list[str | int | float]]:
        temp = []

        if len(diff_tips_between_clade_and_curr_tree) == 0:
            temp.append("monophyletic")
        else:
            temp.append("not_monophyletic")
        temp.append(stats["mean"])
        temp.append(stats["maximum"])
        temp.append(stats["minimum"])
        temp.append(stats["standard_deviation"])
        temp.append(diff_tips_between_clade_and_curr_tree)
        res_arr.append(temp)

        return res_arr

    def print_results(self, res_arr: list[list[str | int | float]]) -> None:
        if self.json_output:
            round_ = round
            sorted_ = sorted
            rows = [
                {
                    "status": res[0],
                    "mean_support": round_(res[1], 4),
                    "max_support": round_(res[2], 4),
                    "min_support": round_(res[3], 4),
                    "stdev_support": round_(res[4], 4),
                    "offending_taxa": (
                        sorted_(res[5]) if len(res) > 5 and res[5] else []
                    ),
                }
                if len(res) > 1
                else {"status": res[0]}
                for res in res_arr
            ]
            print_json(dict(rows=rows, results=rows))
            return

        lines = []
        append = lines.append
        round_ = round
        for res in res_arr:
            try:
                if res[5]:
                    res[5].sort()
                    append(
                        f"{res[0]}\t{round_(res[1], 4)}\t{round_(res[2], 4)}\t{round_(res[3], 4)}\t{round_(res[4], 4)}\t{';'.join(res[5])}"
                    )
                else:
                    append(
                        f"{res[0]}\t{round_(res[1], 4)}\t{round_(res[2], 4)}\t{round_(res[3], 4)}\t{round_(res[4], 4)}"
                    )
            except IndexError:
                append(f"{res[0]}")
        if lines:
            sys.stdout.write("\n".join(lines) + "\n")
