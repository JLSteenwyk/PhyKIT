from __future__ import annotations

import sys

from .base import Tree

def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def calculate_summary_statistics_from_arr(*args, **kwargs):
    from ...helpers.stats_summary import (
        calculate_summary_statistics_from_arr as _calculate_summary_statistics_from_arr,
    )

    return _calculate_summary_statistics_from_arr(*args, **kwargs)


def print_summary_statistics(*args, **kwargs):
    from ...helpers.stats_summary import print_summary_statistics as _print_summary_statistics

    return _print_summary_statistics(*args, **kwargs)


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()


class InternalBranchStats(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"], verbose=parsed["verbose"])
        self.json_output = parsed["json_output"]

    def run(self):
        tree = self.read_tree_file_unmodified()
        stats, lengths_and_names \
            = self.calculate_internal_branch_stats(
                tree, include_names=self.verbose
            )

        if self.json_output:
            if self.verbose:
                rows = [
                    {"length": round(length, 4), "terminals": names}
                    for length, names in lengths_and_names
                ]
                print_json(
                    dict(
                        verbose=True,
                        rows=rows,
                        branches=rows,
                    )
                )
            else:
                print_json(dict(verbose=False, summary=stats))
            return

        if self.verbose:
            try:
                lines = [
                    f"{round(length, 4)} {';'.join(names)}"
                    for length, names in lengths_and_names
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

    def get_internal_branch_lengths(
        self,
        tree: Newick.Tree,
        include_names: bool = True,
    ) -> tuple[
        list[float],
        list[tuple[float, list[str]]]
    ]:
        direct_result = self._get_internal_branch_lengths_direct(
            tree,
            include_names=include_names,
        )
        if direct_result is not None:
            return direct_result

        internal_branch_lengths = []
        lengths_and_names = []
        if not include_names:
            for internal_branch in tree.get_nonterminals():
                if internal_branch.branch_length is not None:
                    internal_branch_lengths.append(internal_branch.branch_length)
            return internal_branch_lengths, lengths_and_names

        clade_term_names = {}
        for clade in tree.find_clades(order="postorder"):
            if clade.is_terminal():
                clade_term_names[id(clade)] = [clade.name]
            else:
                names = []
                for child in clade.clades:
                    names.extend(clade_term_names.get(id(child), []))
                clade_term_names[id(clade)] = names

        for internal_branch in tree.get_nonterminals():
            if internal_branch.branch_length is not None:
                internal_branch_lengths.append(internal_branch.branch_length)
                term_names = clade_term_names.get(id(internal_branch), [])
                lengths_and_names.append(
                    (
                        internal_branch.branch_length, term_names
                    )
                )

        return internal_branch_lengths, lengths_and_names

    @staticmethod
    def _get_internal_branch_lengths_direct(tree, include_names: bool = True):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        internal_branch_lengths = []
        lengths_and_names = []

        if not include_names:
            try:
                stack = [root]
                pop = stack.pop
                append = stack.append
                lengths_append = internal_branch_lengths.append
                while stack:
                    clade = pop()
                    children = clade.clades
                    child_count = len(children)
                    if child_count:
                        if clade.branch_length is not None:
                            lengths_append(clade.branch_length)
                        if child_count == 2:
                            append(children[1])
                            append(children[0])
                        else:
                            for idx in range(child_count - 1, -1, -1):
                                append(children[idx])
            except AttributeError:
                return None
            return internal_branch_lengths, lengths_and_names

        preorder_internal = []
        preorder = []
        try:
            stack = [root]
            pop = stack.pop
            append = stack.append
            preorder_append = preorder.append
            internal_append = preorder_internal.append
            while stack:
                clade = pop()
                preorder_append(clade)
                children = clade.clades
                child_count = len(children)
                if child_count:
                    internal_append(clade)
                    if child_count == 2:
                        append(children[1])
                        append(children[0])
                    else:
                        for idx in range(child_count - 1, -1, -1):
                            append(children[idx])
        except AttributeError:
            return None

        clade_term_names = {}
        try:
            for clade in reversed(preorder):
                children = clade.clades
                child_count = len(children)
                if child_count == 2:
                    names = []
                    names.extend(clade_term_names.get(id(children[0]), []))
                    names.extend(clade_term_names.get(id(children[1]), []))
                    clade_term_names[id(clade)] = names
                elif child_count:
                    names = []
                    for child in children:
                        names.extend(clade_term_names.get(id(child), []))
                    clade_term_names[id(clade)] = names
                else:
                    clade_term_names[id(clade)] = [clade.name]
        except AttributeError:
            return None

        for internal_branch in preorder_internal:
            if internal_branch.branch_length is not None:
                internal_branch_lengths.append(internal_branch.branch_length)
                lengths_and_names.append(
                    (
                        internal_branch.branch_length,
                        clade_term_names.get(id(internal_branch), []),
                    )
                )

        return internal_branch_lengths, lengths_and_names

    @staticmethod
    def _get_internal_branch_lengths_array_direct(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        def iter_lengths():
            stack = [root]
            append = stack.append
            pop = stack.pop
            while stack:
                clade = pop()
                children = clade.clades
                if children:
                    if clade.branch_length is not None:
                        yield clade.branch_length
                    child_count = len(children)
                    if child_count == 2:
                        append(children[1])
                        append(children[0])
                    else:
                        for idx in range(child_count - 1, -1, -1):
                            append(children[idx])

        try:
            return np.fromiter(iter_lengths(), dtype=float)
        except AttributeError:
            return None

    def calculate_internal_branch_stats(
        self,
        tree: Newick.Tree,
        include_names: bool = True,
    ) -> tuple[
        dict[str, float],
        list[tuple[float, list[str]]],
    ]:
        if not include_names:
            internal_branch_lengths_arr = \
                self._get_internal_branch_lengths_array_direct(tree)
            if internal_branch_lengths_arr is not None:
                if internal_branch_lengths_arr.size == 0:
                    print("Calculating internal branch statistics requires a phylogeny with branch lengths.")
                    sys.exit(2)
                stats = calculate_summary_statistics_from_arr(
                    internal_branch_lengths_arr
                )
                return stats, []

        internal_branch_lengths, lengths_and_names = \
            self.get_internal_branch_lengths(tree, include_names=include_names)

        if not internal_branch_lengths:
            print("Calculating internal branch statistics requires a phylogeny with branch lengths.")
            sys.exit(2)

        stats = calculate_summary_statistics_from_arr(internal_branch_lengths)

        return stats, lengths_and_names
