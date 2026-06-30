from __future__ import annotations

from .base import Tree


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class BranchLengthMultiplier(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            factor=parsed["factor"],
            output_file_path=parsed["output_file_path"],
        )
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        if self.factor == 1.0:
            tree = self.read_tree_file_unmodified()
            scaled_count = (
                self.count_branch_lengths(tree) if self.json_output else None
            )
        else:
            tree = self.read_tree_file()
            scaled_count = self.multiply_branch_lengths_by_factor(tree, self.factor)
        self.write_tree_file(tree, self.output_file_path)

        if self.json_output:
            print_json(
                dict(
                    input_tree=self.tree_file_path,
                    factor=self.factor,
                    scaled_branches=scaled_count,
                    output_file=self.output_file_path,
                )
            )

    def process_args(self, args) -> dict[str, str]:
        output_file_path = \
            args.output or f"{args.tree}.factor_{args.factor}.tre"
        return dict(
            tree_file_path=args.tree,
            factor=args.factor,
            output_file_path=output_file_path,
            json_output=getattr(args, "json", False),
        )

    def multiply_branch_lengths_by_factor(
        self,
        tree: Newick.Tree,
        factor: float,
    ) -> int:
        scaled_count = self._multiply_standard_tree_branch_lengths(tree, factor)
        if scaled_count is not None:
            return scaled_count

        scaled_count = 0
        for node in tree.find_clades(order="preorder"):
            if node.branch_length is not None:
                node.branch_length *= factor
                scaled_count += 1
        return scaled_count

    def count_branch_lengths(self, tree: Newick.Tree) -> int:
        count = self._count_standard_tree_branch_lengths(tree)
        if count is not None:
            return count

        branch_count = 0
        for node in tree.find_clades(order="preorder"):
            if node.branch_length is not None:
                branch_count += 1
        return branch_count

    @staticmethod
    def _multiply_standard_tree_branch_lengths(
        tree: Newick.Tree,
        factor: float,
    ):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        stack = [root]
        scaled_count = 0
        try:
            pop = stack.pop
            extend = stack.extend
            while stack:
                node = pop()
                if node.branch_length is not None:
                    node.branch_length *= factor
                    scaled_count += 1
                extend(node.clades)
        except AttributeError:
            return None

        return scaled_count

    @staticmethod
    def _count_standard_tree_branch_lengths(tree: Newick.Tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        stack = [root]
        branch_count = 0
        try:
            pop = stack.pop
            extend = stack.extend
            while stack:
                node = pop()
                if node.branch_length is not None:
                    branch_count += 1
                extend(node.clades)
        except AttributeError:
            return None

        return branch_count
