import sys
from typing import (
    Dict,
    List,
    Tuple,
)

from Bio.Phylo import Newick

from .base import Tree

from ...helpers.stats_summary import calculate_summary_statistics_from_arr, print_summary_statistics
from ...helpers.json_output import print_json


class InternalBranchStats(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"], verbose=parsed["verbose"])
        self.json_output = parsed["json_output"]

    def run(self):
        tree = self.read_tree_file()
        stats, lengths_and_names \
            = self.calculate_internal_branch_stats(tree)

        if self.json_output:
            if self.verbose:
                rows = [
                    dict(length=round(length, 4), terminals=names)
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
                for length, names in lengths_and_names:
                    print(round(length, 4), ";".join(names))
            except BrokenPipeError:
                pass
        else:
            print_summary_statistics(stats)

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            tree_file_path=args.tree,
            verbose=args.verbose,
            json_output=getattr(args, "json", False),
        )

    def get_internal_branch_lengths(
        self,
        tree: Newick.Tree
    ) -> Tuple[
        List[float],
        List[Tuple[float, List[str]]]
    ]:
        internal_branch_lengths = []
        lengths_and_names = []

        # Collect branch lengths and associated names in one pass
        for internal_branch in tree.get_nonterminals():
            if internal_branch.branch_length is not None:
                internal_branch_lengths.append(internal_branch.branch_length)
                term_names = [
                    term.name for term in internal_branch.get_terminals()
                ]
                lengths_and_names.append(
                    (
                        internal_branch.branch_length, term_names
                    )
                )

        return internal_branch_lengths, lengths_and_names

    def calculate_internal_branch_stats(
        self,
        tree: Newick.Tree
    ) -> Tuple[
        Dict[str, float],
        List[Tuple[float, List[str]]],
    ]:
        internal_branch_lengths, lengths_and_names = \
            self.get_internal_branch_lengths(tree)

        if not internal_branch_lengths:
            print("Calculating internal branch statistics requires a phylogeny with branch lengths.")
            sys.exit(2)

        stats = calculate_summary_statistics_from_arr(internal_branch_lengths)

        return stats, lengths_and_names
