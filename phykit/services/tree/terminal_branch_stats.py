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


class TerminalBranchStats(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"], verbose=parsed["verbose"])
        self.json_output = parsed["json_output"]

    def run(self):
        tree = self.read_tree_file_unmodified()
        _, stats, lengths_and_names = \
            self.calculate_terminal_branch_stats(
                tree, include_names=self.verbose
            )

        if self.json_output:
            if self.verbose:
                rows = [
                    dict(length=round(len_and_name[0], 4), taxon=len_and_name[1])
                    for len_and_name in lengths_and_names
                ]
                print_json(
                    dict(
                        verbose=True,
                        rows=rows,
                        tips=rows,
                    )
                )
            else:
                print_json(dict(verbose=False, summary=stats))
            return

        if self.verbose:
            try:
                lines = [
                    f"{round(len_and_name[0], 4)} {len_and_name[1]}"
                    for len_and_name in lengths_and_names
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

    def get_terminal_branch_lengths(
        self,
        tree: Newick.Tree,
        include_names: bool = True,
    ) -> tuple[
        list[float],
        list[list[float | str]],
    ]:
        """
        loop through tree and get all terminal branch lengths
        """
        direct_result = self._get_terminal_branch_lengths_direct(
            tree,
            include_names=include_names,
        )
        if direct_result is not None:
            return direct_result

        terminal_branch_lengths = []
        lengths_and_names = []
        for terminal_branch in tree.get_terminals():
            if terminal_branch.branch_length is not None:
                terminal_branch_lengths.append(terminal_branch.branch_length)
                if include_names:
                    lengths_and_names.append(
                        [terminal_branch.branch_length, terminal_branch.name]
                    )

        return terminal_branch_lengths, lengths_and_names

    @staticmethod
    def _get_terminal_branch_lengths_direct(tree, include_names: bool = True):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        terminal_branch_lengths = []
        lengths_and_names = []
        stack = [root]
        try:
            append_stack = stack.append
            pop = stack.pop
            append_length = terminal_branch_lengths.append
            append_row = lengths_and_names.append
            while stack:
                clade = pop()
                children = clade.clades
                if children:
                    child_count = len(children)
                    if child_count == 2:
                        append_stack(children[1])
                        append_stack(children[0])
                    else:
                        for idx in range(child_count - 1, -1, -1):
                            append_stack(children[idx])
                elif clade.branch_length is not None:
                    branch_length = clade.branch_length
                    append_length(branch_length)
                    if include_names:
                        append_row([branch_length, clade.name])
        except AttributeError:
            return None

        return terminal_branch_lengths, lengths_and_names

    @staticmethod
    def _get_terminal_branch_lengths_array_direct(tree):
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
                    child_count = len(children)
                    if child_count == 2:
                        append(children[1])
                        append(children[0])
                    else:
                        for idx in range(child_count - 1, -1, -1):
                            append(children[idx])
                elif clade.branch_length is not None:
                    yield clade.branch_length

        try:
            return np.fromiter(iter_lengths(), dtype=float)
        except AttributeError:
            return None

    def check_tree_has_branch_lengths(
        self,
        terminal_branch_lengths: list[float],
    ) -> None:
        """
        if tree has no branch lengths, exit
        """
        if len(terminal_branch_lengths) == 0:
            print(
                "Calculating terminal branch statistics requires a phylogeny with branch lengths."
            )
            sys.exit(2)

    def calculate_terminal_branch_stats(
        self,
        tree: Newick.Tree,
        include_names: bool = True,
    ) -> tuple[
        list[float],
        dict[str, float],
        list[list[float | str]],
    ]:
        if not include_names:
            terminal_branch_lengths_arr = \
                self._get_terminal_branch_lengths_array_direct(tree)
            if terminal_branch_lengths_arr is not None:
                if terminal_branch_lengths_arr.size == 0:
                    self.check_tree_has_branch_lengths([])
                stats = calculate_summary_statistics_from_arr(
                    terminal_branch_lengths_arr
                )
                return terminal_branch_lengths_arr.tolist(), stats, []

        terminal_branch_lengths, lengths_and_names = \
            self.get_terminal_branch_lengths(tree, include_names=include_names)

        self.check_tree_has_branch_lengths(terminal_branch_lengths)

        stats = calculate_summary_statistics_from_arr(terminal_branch_lengths)

        return terminal_branch_lengths, stats, lengths_and_names
