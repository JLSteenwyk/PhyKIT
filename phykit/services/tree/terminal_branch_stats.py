import sys
from typing import Dict, List, Tuple, Union

from Bio.Phylo import Newick

from .base import Tree

from ...helpers.stats_summary import (
    calculate_summary_statistics_from_arr,
    print_summary_statistics,
)


class TerminalBranchStats(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        _, stats, lengths_and_names = \
            self.calculate_terminal_branch_stats(tree)

        if self.verbose:
            try:
                for len_and_name in lengths_and_names:
                    print(round(len_and_name[0], 4), len_and_name[1])
            except BrokenPipeError:
                pass
        else:
            print_summary_statistics(stats)

    def process_args(self, args) -> Dict[str, str]:
        return dict(tree_file_path=args.tree, verbose=args.verbose)

    def get_terminal_branch_lengths(
        self,
        tree: Newick.Tree,
    ) -> List[
        Union[
            float,
            List[List[Union[float, str]]],
        ]
    ]:
        """
        loop through tree and get all terminal branch lengths
        """
        terminal_branch_lengths = []
        lengths_and_names = []
        for terminal_branch in tree.get_terminals():
            if terminal_branch.branch_length is not None:
                temp = []
                temp.append(terminal_branch.branch_length)
                terminal_branch_lengths.append(terminal_branch.branch_length)
                temp.append(terminal_branch.name)
                lengths_and_names.append(temp)

        return terminal_branch_lengths, lengths_and_names

    def check_tree_has_branch_lengths(
        self,
        terminal_branch_lengths: List[float],
    ) -> None:
        """
        if tree has no branch lengths, exit
        """
        if len(terminal_branch_lengths) == 0:
            print(
                "Calculating terminal branch statistics requires a phylogeny with branch lengths."
            )
            sys.exit()

    def calculate_terminal_branch_stats(
        self,
        tree: Newick.Tree,
    ) -> Tuple[
        List[float],
        Dict[str, float],
        List[List[Union[float, str]]],
    ]:
        terminal_branch_lengths, lengths_and_names = \
            self.get_terminal_branch_lengths(tree)

        self.check_tree_has_branch_lengths(terminal_branch_lengths)

        stats = calculate_summary_statistics_from_arr(terminal_branch_lengths)

        return terminal_branch_lengths, stats, lengths_and_names
