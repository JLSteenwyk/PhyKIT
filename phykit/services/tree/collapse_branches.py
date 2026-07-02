from __future__ import annotations

from .base import Tree


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class CollapseBranches(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            support=parsed["support"],
            output_file_path=parsed["output_file_path"],
        )
        self.json_output = parsed["json_output"]

    def run(self):
        tree = self.read_tree_file_unmodified()
        copied_tree = False
        initial_nonterminals = None
        has_collapsible_branch = None
        if self.json_output:
            scan_result = self._scan_standard_tree_for_collapse(
                tree,
                self.support,
            )
            if scan_result is not None:
                initial_nonterminals, has_collapsible_branch = scan_result
            else:
                tree = self._fast_copy(tree)
                copied_tree = True
                initial_nonterminals = self.count_internal_nodes(tree)
        else:
            has_collapsible_branch = self._has_collapsible_branch_standard_tree(
                tree,
                self.support,
            )

        if has_collapsible_branch is None:
            if not copied_tree:
                tree = self._fast_copy(tree)
                copied_tree = True
            has_collapsible_branch = self._has_collapsible_branch(tree)
        if has_collapsible_branch is not False:
            if not copied_tree:
                tree = self._fast_copy(tree)
            tree.collapse_all(
                lambda c: c.confidence and c.confidence < self.support
            )
        self.write_tree_file(tree, self.output_file_path)

        if self.json_output:
            final_nonterminals = (
                initial_nonterminals
                if has_collapsible_branch is False
                else self.count_internal_nodes(tree)
            )
            print_json(
                dict(
                    input_tree=self.tree_file_path,
                    support_threshold=self.support,
                    collapsed_branches=max(0, initial_nonterminals - final_nonterminals),
                    output_file=self.output_file_path,
                )
            )

    def process_args(self, args) -> dict[str, str]:
        output_file_path = \
            args.output or f"{args.tree}.collapsed_{args.support}.tre"
        return dict(
            tree_file_path=args.tree,
            support=args.support,
            output_file_path=output_file_path,
            json_output=getattr(args, "json", False),
        )

    def count_internal_nodes(self, tree) -> int:
        count = self._count_standard_tree_internal_nodes(tree)
        if count is not None:
            return count
        return len(tree.get_nonterminals())

    def _has_collapsible_branch(self, tree):
        result = self._has_collapsible_branch_standard_tree(tree, self.support)
        if result is not None:
            return result
        return None

    @staticmethod
    def _has_collapsible_branch_standard_tree(tree, support):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        stack = [root]
        try:
            pop = stack.pop
            extend = stack.extend
            while stack:
                node = pop()
                confidence = node.confidence
                if confidence and confidence < support:
                    return True
                children = node.clades
                if children:
                    extend(children)
        except AttributeError:
            return None
        return False

    @staticmethod
    def _scan_standard_tree_for_collapse(tree, support):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        internal_count = 0
        has_collapsible_branch = False
        stack = [root]
        try:
            while stack:
                node = stack.pop()
                children = node.clades
                if children:
                    internal_count += 1
                confidence = node.confidence
                if confidence and confidence < support:
                    has_collapsible_branch = True
                if children:
                    stack.extend(children)
        except AttributeError:
            return None
        return internal_count, has_collapsible_branch

    @staticmethod
    def _count_standard_tree_internal_nodes(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        count = 0
        stack = [root]
        try:
            while stack:
                node = stack.pop()
                children = node.clades
                if children:
                    count += 1
                    stack.extend(children)
        except AttributeError:
            return None
        return count
