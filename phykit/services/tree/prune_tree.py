from __future__ import annotations

import re

from .base import Tree


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def read_single_column_file_to_list(*args, **kwargs):
    from ...helpers.files import (
        read_single_column_file_to_list as _read_single_column_file_to_list,
    )

    return _read_single_column_file_to_list(*args, **kwargs)


_BRANCH_LABEL_RE = re.compile(r"\{[^{}]*\}")


def _strip_branch_label(name: str) -> str:
    """Remove HyPhy/aBSREL-style {…} branch labels from a tip name."""
    return _BRANCH_LABEL_RE.sub("", name) if name and "{" in name else name


class PruneTree(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            list_of_taxa=parsed["list_of_taxa"],
            output_file_path=parsed["output_file_path"],
            keep=parsed["keep"],
        )
        self.json_output = parsed["json_output"]
        self.ignore_branch_labels = parsed["ignore_branch_labels"]

    def run(self) -> None:
        tree = self.read_tree_file_unmodified()

        taxa = read_single_column_file_to_list(self.list_of_taxa)

        if self.ignore_branch_labels:
            taxa_set = set(taxa)
            tips_in_tree = self.get_tip_names_from_tree(tree)
            if self.keep:
                taxa = [
                    tip for tip in tips_in_tree
                    if _strip_branch_label(tip) not in taxa_set
                ]
            else:
                taxa = [
                    tip for tip in tips_in_tree
                    if _strip_branch_label(tip) in taxa_set
                ]
        elif self.keep:
            taxa_set = set(taxa)
            tips_in_tree = self.get_tip_names_from_tree(tree)
            taxa = [x for x in tips_in_tree if x not in taxa_set]

        if taxa:
            tree = self._fast_copy(tree)
            tree = self.prune_tree_using_taxa_list(tree, taxa)

        self.write_tree_file(tree, self.output_file_path)

        if self.json_output:
            remaining_tips = self.calculate_terminal_count_fast(tree)
            if remaining_tips is None:
                remaining_tips = tree.count_terminals()
            print_json(
                dict(
                    input_tree=self.tree_file_path,
                    input_taxa_file=self.list_of_taxa,
                    keep_input_taxa=self.keep,
                    ignore_branch_labels=self.ignore_branch_labels,
                    taxa_pruned=sorted(taxa),
                    pruned_count=len(taxa),
                    remaining_tips=remaining_tips,
                    output_file=self.output_file_path,
                )
            )

    def process_args(self, args) -> dict[str, str]:
        tree_file_path = args.tree
        output_file_path = \
            f"{args.output}" if args.output else f"{tree_file_path}.pruned"

        keep = True if args.keep is None else args.keep

        return dict(
            tree_file_path=tree_file_path,
            list_of_taxa=args.list_of_taxa,
            output_file_path=output_file_path,
            keep=keep,
            ignore_branch_labels=getattr(args, "ignore_branch_labels", False),
            json_output=getattr(args, "json", False),
        )
