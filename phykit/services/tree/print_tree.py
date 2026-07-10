from __future__ import annotations

from .base import Tree


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class _LazyPhylo:
    def draw_ascii(self, *args, **kwargs):
        from Bio import Phylo as _Phylo

        self.draw_ascii = _Phylo.draw_ascii
        return self.draw_ascii(*args, **kwargs)


Phylo = _LazyPhylo()


class PrintTree(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            remove=parsed["remove"],
        )
        self.json_output = parsed["json_output"]

    def run(self):
        tree = self.read_tree_file() if self.remove else self.read_tree_file_unmodified()

        if self.remove:
            self.remove_branch_lengths(tree)

        if self.json_output:
            print_json(
                dict(
                    remove_branch_lengths=self.remove,
                    tree_newick=tree.format("newick").strip(),
                )
            )
            return

        try:
            Phylo.draw_ascii(tree)
        except BrokenPipeError:
            pass

    def process_args(self, args) -> dict[str, str]:
        return dict(
            tree_file_path=args.tree,
            remove=args.remove,
            json_output=getattr(args, "json", False),
        )

    def remove_branch_lengths(self, tree) -> None:
        if self._remove_standard_tree_branch_lengths(tree):
            return

        for node in tree.find_clades(order="preorder"):
            node.branch_length = None

    @staticmethod
    def _remove_standard_tree_branch_lengths(tree) -> bool:
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return False

        stack = [root]
        try:
            pop = stack.pop
            extend = stack.extend
            while stack:
                node = pop()
                node.branch_length = None
                extend(node.clades)
        except AttributeError:
            return False
        return True
