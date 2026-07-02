"""Utilities for deriving paths through Bio.Phylo trees."""

from __future__ import annotations


def build_object_parent_map(tree) -> dict[object, object]:
    """Build a mapping from child clade object to parent clade object."""
    direct_parent_map = _build_object_parent_map_direct(tree)
    if direct_parent_map is not None:
        return direct_parent_map

    parent_map = {}
    for clade in tree.find_clades(order="preorder"):
        for child in clade.clades:
            parent_map[child] = clade
    return parent_map


def _build_object_parent_map_direct(tree):
    try:
        root = tree.root
        root.clades
    except AttributeError:
        return None

    parent_map = {}
    stack = [root]
    try:
        pop = stack.pop
        extend = stack.extend
        while stack:
            clade = pop()
            children = clade.clades
            for child in children:
                parent_map[child] = clade
            if children:
                extend(children)
    except AttributeError:
        return None
    return parent_map


def path_from_root(
    node,
    root,
    parent_map: dict[object, object],
) -> list[object] | None:
    """Return the path from root to node, or None when parent links are incomplete."""
    path = [node]
    current = node
    remaining = len(parent_map) + 1

    while current is not root:
        if remaining <= 0:
            return None
        parent = parent_map.get(current)
        if parent is None:
            return None
        current = parent
        path.append(current)
        remaining -= 1

    path.reverse()
    return path


def build_root_path_map(tree, root=None) -> dict[object, list[object]]:
    """Build root-to-node paths for every clade in one preorder traversal."""
    if root is None:
        root = tree.root
    elif root is not tree.root:
        return _build_root_path_map_legacy(tree, root)

    direct_paths = _build_root_path_map_direct(root)
    if direct_paths is not None:
        return direct_paths

    return _build_root_path_map_legacy(tree, root)


def _build_root_path_map_legacy(tree, root) -> dict[object, list[object]]:
    paths = {root: [root]}
    for clade in tree.find_clades(order="preorder"):
        parent_path = paths.get(clade)
        if parent_path is None:
            return {}
        for child in clade.clades:
            paths[child] = parent_path + [child]
    return paths


def _build_root_path_map_direct(root):
    try:
        root.clades
    except AttributeError:
        return None

    paths = {root: [root]}
    stack = [(root, paths[root])]
    try:
        while stack:
            clade, path = stack.pop()
            children = clade.clades
            for child in reversed(children):
                child_path = path + [child]
                paths[child] = child_path
                stack.append((child, child_path))
    except AttributeError:
        return None
    return paths
