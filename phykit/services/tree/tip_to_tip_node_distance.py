from __future__ import annotations

import sys

from .base import Tree


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class _LazyTreeMixin:
    def find_any(self, *args, **kwargs):
        from Bio.Phylo.BaseTree import TreeMixin as _TreeMixin

        return _TreeMixin.find_any(*args, **kwargs)

    def trace(self, *args, **kwargs):
        from Bio.Phylo.BaseTree import TreeMixin as _TreeMixin

        return _TreeMixin.trace(*args, **kwargs)


TreeMixin = _LazyTreeMixin()


class TipToTipNodeDistance(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            tip_1=parsed["tip_1"],
            tip_2=parsed["tip_2"],
        )
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        tree_zero = self.read_tree_file_unmodified()

        node_distance = self.calculate_tip_to_tip_node_distance(
            tree_zero, self.tip_1, self.tip_2
        )
        if self.json_output:
            print_json(
                dict(
                    taxon_a=self.tip_1,
                    taxon_b=self.tip_2,
                    tip_to_tip_node_distance=node_distance,
                )
            )
            return
        print(node_distance)

    def check_leaves(
        self,
        tree_zero: Newick.Tree,
        tip_1: str,
        tip_2: str,
    ) -> None:
        leaf1 = TreeMixin.find_any(tree_zero, tip_1)
        if not bool(leaf1):
            print(tip_1, "not on tree\nExiting...")
            sys.exit(2)
        leaf2 = TreeMixin.find_any(tree_zero, tip_2)
        if not bool(leaf2):
            print(tip_2, "not on tree\nExiting...")
            sys.exit(2)

    def calculate_tip_to_tip_node_distance(
        self,
        tree_zero: Newick.Tree,
        tip_1: str,
        tip_2: str,
    ) -> int:
        node_distance = self._calculate_tip_to_tip_node_distance_fast(
            tree_zero, tip_1, tip_2
        )
        if node_distance is not None:
            return node_distance

        self.check_leaves(tree_zero, tip_1, tip_2)
        return len(TreeMixin.trace(tree_zero, tip_1, tip_2))

    @staticmethod
    def _calculate_tip_to_tip_node_distance_fast(
        tree_zero: Newick.Tree,
        tip_1: str,
        tip_2: str,
    ):
        try:
            root = tree_zero.root
            if tip_1 == tip_2:
                return TipToTipNodeDistance._calculate_same_tip_node_distance_fast(
                    root,
                    tip_1,
                )

            parent_map = {}
            found = {}
            targets = {tip_1, tip_2}
            stack = [root]
            pop = stack.pop
            append = stack.append

            while stack and len(found) < len(targets):
                clade = pop()
                children = clade.clades
                if children:
                    child_count = len(children)
                    if child_count == 2:
                        left, right = children
                        parent_map[right] = clade
                        append(right)
                        parent_map[left] = clade
                        append(left)
                    else:
                        for index in range(child_count - 1, -1, -1):
                            child = children[index]
                            parent_map[child] = clade
                            append(child)
                elif clade.name in targets:
                    found.setdefault(clade.name, clade)

            leaf1 = found.get(tip_1)
            if leaf1 is None:
                print(tip_1, "not on tree\nExiting...")
                sys.exit(2)

            leaf2 = found.get(tip_2)
            if leaf2 is None:
                print(tip_2, "not on tree\nExiting...")
                sys.exit(2)

            root = tree_zero.root
            ancestors = {}
            current = leaf1
            distance = 0
            while True:
                ancestors[current] = distance
                if current is root:
                    break
                current = parent_map[current]
                distance += 1

            current = leaf2
            distance = 0
            while current not in ancestors:
                current = parent_map[current]
                distance += 1

            return ancestors[current] + distance
        except (AttributeError, KeyError, TypeError):
            return None

    @staticmethod
    def _calculate_same_tip_node_distance_fast(root, tip_name: str):
        stack = [root]
        pop = stack.pop
        append = stack.append
        try:
            while stack:
                clade = pop()
                children = clade.clades
                if children:
                    child_count = len(children)
                    if child_count == 2:
                        append(children[1])
                        append(children[0])
                    else:
                        for index in range(child_count - 1, -1, -1):
                            append(children[index])
                elif clade.name == tip_name:
                    return 0
        except AttributeError:
            return None

        print(tip_name, "not on tree\nExiting...")
        sys.exit(2)

    @staticmethod
    def _path_to_root(node, root, parent_map):
        path = [node]
        current = node
        while current is not root:
            current = parent_map[current]
            path.append(current)
        path.reverse()
        return path

    def process_args(self, args) -> dict[str, str]:
        return dict(
            tree_file_path=args.tree_zero,
            tip_1=args.tip_1,
            tip_2=args.tip_2,
            json_output=getattr(args, "json", False),
        )
