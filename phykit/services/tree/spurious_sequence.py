import statistics as stat
from typing import Dict, List, Tuple

from Bio.Phylo import Newick

from .base import Tree


class SpuriousSequence(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        tree = self.read_tree_file()
        name_and_branch_len, threshold, median = \
            self.identify_spurious_sequence(
                tree, self.factor
            )

        counter = 0
        for name, length in name_and_branch_len.items():
            if length >= threshold:
                try:
                    print(
                        f"{name}\t{round(length, 4)}\t{round(threshold, 4)}\t{round(median, 4)}"
                    )
                except BrokenPipeError:
                    pass
                counter += 1

        if counter == 0:
            print("None")

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            tree_file_path=args.tree,
            factor=args.factor or 20
        )

    def identify_spurious_sequence(
        self,
        tree: Newick.Tree,
        factor: float,
    ) -> Tuple[
        Dict[str, float],
        float,
        float
    ]:
        branch_lengths, name_and_branch_len = \
            self.get_branch_lengths_and_their_names(tree)

        median = stat.median(branch_lengths)

        threshold = median * factor

        return name_and_branch_len, threshold, median

    def get_branch_lengths_and_their_names(
        self,
        tree: Newick.Tree,
    ) -> Tuple[
        List[float],
        Dict[str, float],
    ]:
        branch_lengths = list()
        name_and_branch_len = dict()

        # collect terminal branch lengths
        for terminal in tree.get_terminals():
            branch_lengths.append(terminal.branch_length)
            name_and_branch_len[terminal.name] = terminal.branch_length

        # collect internal branch lengths
        for internal in tree.get_nonterminals():
            branch_lengths.append(terminal.branch_length)

        return branch_lengths, name_and_branch_len
