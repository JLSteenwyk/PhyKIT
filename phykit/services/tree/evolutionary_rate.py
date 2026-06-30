from __future__ import annotations

from .base import Tree


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class EvolutionaryRate(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        tree = self.read_tree_file_unmodified()
        total_tree_length, num_terminals = (
            self.calculate_total_branch_length_and_terminal_count_fast(tree)
        )
        evo_rate = round(total_tree_length / num_terminals, 4)
        if self.json_output:
            print_json(dict(evolutionary_rate=evo_rate))
            return
        print(evo_rate)

    def process_args(self, args) -> dict[str, str]:
        return dict(tree_file_path=args.tree, json_output=getattr(args, "json", False))
