from __future__ import annotations

from .base import Tree


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class TotalTreeLength(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        summary = self._get_simple_newick_summary(
            self.tree_file_path,
            "tree_file_path",
        )
        if summary is None:
            tree = self.read_tree_file_unmodified()
            total_tree_length = round(self.calculate_total_tree_length(tree), 4)
        else:
            _, total_len, _ = summary
            total_tree_length = round(total_len, 4)
        if self.json_output:
            print_json(dict(total_tree_length=total_tree_length))
            return
        print(total_tree_length)

    def process_args(self, args) -> dict[str, str]:
        return dict(tree_file_path=args.tree, json_output=getattr(args, "json", False))

    def calculate_total_tree_length(
        self,
        tree: Newick.Tree
    ) -> int | float:
        total_len = self.calculate_total_branch_length_fast(tree)

        if isinstance(total_len, (int, float)):
            return total_len
