from __future__ import annotations

import sys

from .base import Tree


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class RenameTreeTips(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            idmap=parsed["idmap"],
            output_file_path=parsed["output_file_path"],
        )
        self.json_output = parsed["json_output"]

    def run(self):
        idmap = self.read_id_map()

        tree = self.read_tree_file_unmodified()
        if self.json_output:
            should_rename = self.count_matching_tip_names(tree, idmap) > 0
        else:
            should_rename = self.has_matching_tip_name(tree, idmap)

        renamed_count = 0
        if should_rename:
            tree = self._fast_copy(tree)
            tree, renamed_count = self.replace_tip_names(tree, idmap)

        self.write_tree_file(tree, self.output_file_path)

        if self.json_output:
            print_json(
                dict(
                    input_tree=self.tree_file_path,
                    idmap=self.idmap,
                    output_file=self.output_file_path,
                    renamed_tips=renamed_count,
                )
            )

    def process_args(self, args) -> dict[str, str]:
        tree_file_path = args.tree

        output_file_path = \
            f"{args.output}" if args.output else f"{tree_file_path}.renamed"

        return dict(
            tree_file_path=tree_file_path,
            idmap=args.idmap,
            output_file_path=output_file_path,
            json_output=getattr(args, "json", False),
        )

    def read_id_map(self) -> dict[str, str]:
        idmap = dict()
        try:
            with open(self.idmap, "rb") as identifiers:
                for line in identifiers:
                    (key, val) = line.split()
                    idmap[key.decode()] = val.decode()
        except FileNotFoundError:
            try:
                print(f"{self.idmap} corresponds to no such file.")
                print("Please check file name and pathing")
                sys.exit(2)
            except BrokenPipeError:
                pass

        return idmap

    def replace_tip_names(
        self,
        tree: Tree,
        idmap: dict[str, str]
    ) -> tuple[Newick.Tree, int]:
        if not idmap:
            return tree, 0

        result = self._replace_standard_tree_tip_names(tree, idmap)
        if result is not None:
            return result

        renamed_count = 0
        for term in tree.get_terminals():
            name = term.name
            if name in idmap:
                term.name = idmap[name]
                renamed_count += 1

        return tree, renamed_count

    def count_matching_tip_names(self, tree: Tree, idmap: dict[str, str]) -> int:
        if not idmap:
            return 0

        count = self._count_standard_tree_matching_tip_names(tree, idmap)
        if count is not None:
            return count

        renamed_count = 0
        for term in tree.get_terminals():
            if term.name in idmap:
                renamed_count += 1
        return renamed_count

    def has_matching_tip_name(self, tree: Tree, idmap: dict[str, str]) -> bool:
        if not idmap:
            return False

        result = self._has_standard_tree_matching_tip_name(tree, idmap)
        if result is not None:
            return result

        return any(term.name in idmap for term in tree.get_terminals())

    @staticmethod
    def _replace_standard_tree_tip_names(
        tree: Tree,
        idmap: dict[str, str],
    ):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        stack = [root]
        renamed_count = 0
        try:
            pop = stack.pop
            extend = stack.extend
            contains = idmap.__contains__
            get = idmap.__getitem__
            while stack:
                node = pop()
                children = node.clades
                if children:
                    extend(children)
                else:
                    name = node.name
                    if contains(name):
                        node.name = get(name)
                        renamed_count += 1
        except AttributeError:
            return None

        return tree, renamed_count

    @staticmethod
    def _count_standard_tree_matching_tip_names(
        tree: Tree,
        idmap: dict[str, str],
    ):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        renamed_count = 0
        stack = [root]
        pop = stack.pop
        extend = stack.extend
        contains = idmap.__contains__
        try:
            while stack:
                node = pop()
                children = node.clades
                if children:
                    extend(children)
                elif contains(node.name):
                    renamed_count += 1
        except AttributeError:
            return None

        return renamed_count

    @staticmethod
    def _has_standard_tree_matching_tip_name(
        tree: Tree,
        idmap: dict[str, str],
    ):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        stack = [root]
        pop = stack.pop
        extend = stack.extend
        contains = idmap.__contains__
        try:
            while stack:
                node = pop()
                children = node.clades
                if children:
                    extend(children)
                elif contains(node.name):
                    return True
        except AttributeError:
            return None

        return False
