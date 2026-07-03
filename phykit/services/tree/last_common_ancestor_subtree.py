from __future__ import annotations

import sys

from .base import Tree


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def read_single_column_file_to_list(*args, **kwargs):
    from ...helpers.files import (
        read_single_column_file_to_list as _read_single_column_file_to_list,
    )

    return _read_single_column_file_to_list(*args, **kwargs)


def _find_parent_depth_lca(targets, parent_by_clade, depth_by_clade):
    lca = targets[0]
    lca_depth = depth_by_clade[lca]
    for idx in range(1, len(targets)):
        current = targets[idx]
        current_depth = depth_by_clade[current]

        while lca_depth > current_depth:
            lca = parent_by_clade[lca]
            lca_depth -= 1
        while current_depth > lca_depth:
            current = parent_by_clade[current]
            current_depth -= 1
        while lca is not current:
            lca = parent_by_clade[lca]
            current = parent_by_clade[current]
            lca_depth -= 1
        if lca_depth == 0:
            return lca

    return lca


class LastCommonAncestorSubtree(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            output_file_path=parsed["output_file_path"],
            list_of_taxa=parsed["list_of_taxa"],
        )
        self.json_output = parsed["json_output"]

    def run(self):
        tree = self.read_tree_file_unmodified()
        try:
            taxa = read_single_column_file_to_list(self.list_of_taxa)
        except FileNotFoundError:
            print("Taxa list file is not found. Please check pathing.")
            sys.exit(2)
        subtree = self._find_lca_subtree(tree, taxa)

        self.write_tree_file(subtree, self.output_file_path)

        if self.json_output:
            subtree_tip_count = self.calculate_terminal_count_fast(subtree)
            if subtree_tip_count is None:
                subtree_tip_count = subtree.count_terminals()
            print_json(
                dict(
                    input_tree=self.tree_file_path,
                    input_taxa_file=self.list_of_taxa,
                    taxa_count=len(taxa),
                    subtree_tips=subtree_tip_count,
                    output_file=self.output_file_path,
                )
            )

    def process_args(self, args) -> dict[str, str]:
        tree_file_path = args.tree
        output_file_path = args.output or f"{tree_file_path}.subtree.tre"

        return dict(
            tree_file_path=tree_file_path,
            output_file_path=output_file_path,
            list_of_taxa=args.list_of_taxa,
            json_output=getattr(args, "json", False),
        )

    @staticmethod
    def _find_lca_subtree(tree, taxa):
        """Find the MRCA clade for taxa, using parent links when possible."""
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return tree.common_ancestor(taxa)

        selected_names = set(taxa)
        selected_count = len(selected_names)
        if selected_count == 1:
            single_lca = LastCommonAncestorSubtree._find_single_terminal_subtree(
                tree,
                taxa,
                next(iter(selected_names)),
            )
            if single_lca is not None:
                return single_lca

        terminal_by_name = {}
        parent_by_clade = {root: None}
        depth_by_clade = {root: 0}
        stack = [root]
        pop = stack.pop
        append = stack.append
        try:
            while stack:
                clade = pop()
                children = clade.clades
                if children:
                    child_depth = depth_by_clade[clade] + 1
                    n_children = len(children)
                    if n_children == 2:
                        left, right = children
                        parent_by_clade[right] = clade
                        depth_by_clade[right] = child_depth
                        append(right)
                        parent_by_clade[left] = clade
                        depth_by_clade[left] = child_depth
                        append(left)
                    else:
                        for idx in range(n_children - 1, -1, -1):
                            child = children[idx]
                            parent_by_clade[child] = clade
                            depth_by_clade[child] = child_depth
                            append(child)
                else:
                    try:
                        if clade.name in selected_names:
                            terminal_by_name[clade.name] = clade
                            if (
                                selected_count
                                and len(terminal_by_name) == selected_count
                            ):
                                stack.clear()
                                break
                    except TypeError:
                        pass
            targets = [terminal_by_name[taxon] for taxon in taxa]
        except (AttributeError, KeyError, TypeError):
            return tree.common_ancestor(taxa)

        if not targets:
            return tree.common_ancestor(taxa)

        return _find_parent_depth_lca(targets, parent_by_clade, depth_by_clade)

    @staticmethod
    def _find_single_terminal_subtree(tree, taxa, selected_name):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return tree.common_ancestor(taxa)

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
                        for idx in range(child_count - 1, -1, -1):
                            append(children[idx])
                elif clade.name == selected_name:
                    return clade
        except AttributeError:
            return tree.common_ancestor(taxa)

        return tree.common_ancestor(taxa)
