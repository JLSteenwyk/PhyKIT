from __future__ import annotations

import sys

from .base import Tree


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class TipLabels(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        simple_tips = self._get_simple_newick_tip_names(
            self.tree_file_path,
            "tree_file_path",
        )
        if simple_tips is None:
            tree = self.read_tree_file_unmodified()
            tips = self.get_tip_names_from_tree(tree)
        else:
            tips = list(simple_tips)
        if self.json_output:
            rows = [{"taxon": tip} for tip in tips]
            print_json(dict(rows=rows, tips=tips))
            return
        try:
            sys.stdout.write("\n".join(tips) + "\n")
        except BrokenPipeError:
            pass

    def process_args(self, args) -> dict[str, str]:
        return dict(tree_file_path=args.tree, json_output=getattr(args, "json", False))
