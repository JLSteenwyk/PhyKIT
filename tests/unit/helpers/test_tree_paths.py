from io import StringIO
import subprocess
import sys

from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin
from Bio.Phylo.Newick import Clade, Tree as NewickTree

from phykit.helpers.tree_paths import (
    build_object_parent_map,
    build_root_path_map,
    path_from_root,
)


def test_module_import_does_not_import_typing():
    code = """
import sys
import phykit.helpers.tree_paths
assert "typing" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_build_object_parent_map_uses_direct_traversal(monkeypatch):
    tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")

    def fail_find_clades(*_args, **_kwargs):
        raise AssertionError("generic tree traversal should not be used")

    monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

    parent_map = build_object_parent_map(tree)

    left = tree.root.clades[0]
    right = tree.root.clades[1]
    assert parent_map[left] is tree.root
    assert parent_map[right] is tree.root
    assert parent_map[left.clades[0]] is left
    assert parent_map[left.clades[1]] is left
    assert parent_map[right.clades[0]] is right
    assert parent_map[right.clades[1]] is right

    assert path_from_root(right.clades[1], tree.root, parent_map) == [
        tree.root,
        right,
        right.clades[1],
    ]


def test_build_object_parent_map_handles_mixed_child_counts(monkeypatch):
    tree = NewickTree(
        root=Clade(
            clades=[
                Clade(name="A"),
                Clade(clades=[Clade(name="B"), Clade(name="C")]),
                Clade(clades=[Clade(name="D"), Clade(name="E"), Clade(name="F")]),
            ],
        )
    )

    def fail_find_clades(*_args, **_kwargs):
        raise AssertionError("generic tree traversal should not be used")

    monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

    parent_map = build_object_parent_map(tree)

    binary = tree.root.clades[1]
    trifurcation = tree.root.clades[2]
    assert parent_map[tree.root.clades[0]] is tree.root
    assert parent_map[binary] is tree.root
    assert parent_map[trifurcation] is tree.root
    assert parent_map[binary.clades[0]] is binary
    assert parent_map[binary.clades[1]] is binary
    assert parent_map[trifurcation.clades[0]] is trifurcation
    assert parent_map[trifurcation.clades[1]] is trifurcation
    assert parent_map[trifurcation.clades[2]] is trifurcation


def test_build_root_path_map_uses_direct_traversal(monkeypatch):
    tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")

    def fail_find_clades(*_args, **_kwargs):
        raise AssertionError("generic tree traversal should not be used")

    monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

    paths = build_root_path_map(tree)

    left = tree.root.clades[0]
    right = tree.root.clades[1]
    assert paths[tree.root] == [tree.root]
    assert paths[left] == [tree.root, left]
    assert paths[right.clades[1]] == [tree.root, right, right.clades[1]]


def test_build_root_path_map_binary_children_do_not_call_reversed():
    class NoReversedList(list):
        def __reversed__(self):
            raise AssertionError("binary root-path traversal should push explicitly")

    def wrap_binary_children(clade):
        if len(clade.clades) == 2:
            clade.clades = NoReversedList(clade.clades)
        for child in clade.clades:
            wrap_binary_children(child)

    tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
    wrap_binary_children(tree.root)

    paths = build_root_path_map(tree)

    left = tree.root.clades[0]
    right = tree.root.clades[1]
    assert paths[left.clades[0]] == [tree.root, left, left.clades[0]]
    assert paths[right.clades[1]] == [tree.root, right, right.clades[1]]


def test_path_from_root_returns_none_for_incomplete_parent_map():
    root = Clade(name="root")
    parent = Clade(name="parent")
    child = Clade(name="child")

    assert path_from_root(child, root, {child: parent}) is None


def test_path_from_root_returns_none_for_parent_cycle():
    root = Clade(name="root")
    parent = Clade(name="parent")
    child = Clade(name="child")

    assert path_from_root(child, root, {child: parent, parent: child}) is None
