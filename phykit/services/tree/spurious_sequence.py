import statistics as stat
from typing import Dict, List, Tuple

from Bio.Phylo import Newick

from .base import Tree
from ...helpers.json_output import print_json


class SpuriousSequence(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"], factor=parsed["factor"])
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        tree = self.read_tree_file()
        name_and_branch_len, threshold, median = \
            self.identify_spurious_sequence(
                tree, self.factor
            )

        spurious_rows = []
        counter = 0
        for name, length in name_and_branch_len.items():
            if length >= threshold:
                spurious_rows.append(
                    dict(
                        taxon=name,
                        branch_length=round(length, 4),
                        threshold=round(threshold, 4),
                        median=round(median, 4),
                    )
                )
                if not self.json_output:
                    try:
                        print(
                            f"{name}\t{round(length, 4)}\t{round(threshold, 4)}\t{round(median, 4)}"
                        )
                    except BrokenPipeError:
                        pass
                counter += 1

        if self.json_output:
            print_json(dict(rows=spurious_rows, spurious_sequences=spurious_rows))
            return

        if counter == 0:
            print("None")

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            tree_file_path=args.tree,
            factor=args.factor or 20,
            json_output=getattr(args, "json", False),
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
        branch_lengths = []
        name_and_branch_len = {}

        # collect terminal branch lengths only for spurious sequence detection
        # (internal branches are not considered for spurious sequence detection)
        for terminal in tree.get_terminals():
            if terminal.branch_length is not None:
                branch_lengths.append(terminal.branch_length)
                name_and_branch_len[terminal.name] = terminal.branch_length

        return branch_lengths, name_and_branch_len
