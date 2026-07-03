from argparse import Namespace
from io import StringIO
import subprocess
import sys

import pytest
from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin

import phykit.services.tree.last_common_ancestor_subtree as module
from phykit.services.tree.last_common_ancestor_subtree import LastCommonAncestorSubtree


def test_module_import_does_not_import_json_helpers():
    code = """
import sys
import phykit.services.tree.last_common_ancestor_subtree as module

assert callable(module.print_json)
assert callable(module.read_single_column_file_to_list)
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.files" not in sys.modules
assert "typing" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


@pytest.fixture
def args():
    return Namespace(
        tree="/some/path/to/file.tre",
        list_of_taxa="/some/path/to/taxa.txt",
        output=None,
    )


class _Subtree:
    def count_terminals(self):
        return 3


class _Tree:
    def __init__(self):
        self.common_ancestor_calls = []

    def common_ancestor(self, taxa):
        self.common_ancestor_calls.append(taxa)
        assert taxa == ["a", "b"]
        return _Subtree()


class _Clade:
    def __init__(self, name=None, clades=None):
        self.name = name
        self.clades = clades or []


class _ExplodingClade:
    name = "unused"

    @property
    def clades(self):
        raise AssertionError("unneeded subtree should not be traversed")


class _DirectTree:
    def __init__(self, root):
        self.root = root


def test_find_parent_depth_lca_does_not_slice_targets():
    class NoSliceList(list):
        def __getitem__(self, key):
            if isinstance(key, slice):
                raise AssertionError("parent-depth LCA target scan should not slice")
            return super().__getitem__(key)

    root = object()
    internal = object()
    left = object()
    right = object()
    parent_by_clade = {
        root: None,
        internal: root,
        left: internal,
        right: root,
    }
    depth_by_clade = {
        root: 0,
        internal: 1,
        left: 2,
        right: 1,
    }

    assert (
        module._find_parent_depth_lca(
            NoSliceList([left, right]),
            parent_by_clade,
            depth_by_clade,
        )
        is root
    )


def test_find_parent_depth_lca_returns_when_lca_reaches_root():
    root = object()
    left = object()
    right = object()
    unindexed_target = object()
    parent_by_clade = {
        root: None,
        left: root,
        right: root,
    }
    depth_by_clade = {
        root: 0,
        left: 1,
        right: 1,
    }

    assert (
        module._find_parent_depth_lca(
            [left, right, unindexed_target],
            parent_by_clade,
            depth_by_clade,
        )
        is root
    )


class TestLastCommonAncestorSubtree:
    def test_init_sets_expected_attrs(self, args):
        service = LastCommonAncestorSubtree(args)
        assert service.tree_file_path == args.tree
        assert service.list_of_taxa == args.list_of_taxa
        assert service.output_file_path == f"{args.tree}.subtree.tre"
        assert service.json_output is False

    def test_process_args_honors_output_and_json(self):
        args = Namespace(
            tree="/some/path/to/file.tre",
            list_of_taxa="/some/path/to/taxa.txt",
            output="/tmp/subtree.tre",
            json=True,
        )
        service = LastCommonAncestorSubtree(args)
        assert service.output_file_path == "/tmp/subtree.tre"
        assert service.json_output is True

    def test_run_writes_subtree_and_json(self, mocker):
        args = Namespace(
            tree="/some/path/to/file.tre",
            list_of_taxa="/some/path/to/taxa.txt",
            output="/tmp/subtree.tre",
            json=True,
        )
        service = LastCommonAncestorSubtree(args)
        tree = _Tree()
        mocker.patch.object(
            LastCommonAncestorSubtree,
            "read_tree_file_unmodified",
            return_value=tree,
        )
        mocker.patch(
            "phykit.services.tree.last_common_ancestor_subtree.read_single_column_file_to_list",
            return_value=["a", "b"],
        )
        mocked_write = mocker.patch.object(
            LastCommonAncestorSubtree, "write_tree_file"
        )
        mocked_json = mocker.patch(
            "phykit.services.tree.last_common_ancestor_subtree.print_json"
        )

        service.run()

        assert tree.common_ancestor_calls == [["a", "b"]]
        mocked_write.assert_called_once()
        payload = mocked_json.call_args.args[0]
        assert payload["taxa_count"] == 2
        assert payload["subtree_tips"] == 3
        assert payload["output_file"] == "/tmp/subtree.tre"

    def test_run_uses_unmodified_tree_read(self, mocker):
        args = Namespace(
            tree="/some/path/to/file.tre",
            list_of_taxa="/some/path/to/taxa.txt",
            output="/tmp/subtree.tre",
            json=False,
        )
        service = LastCommonAncestorSubtree(args)
        tree = _Tree()
        mocked_read = mocker.patch.object(
            LastCommonAncestorSubtree,
            "read_tree_file_unmodified",
            return_value=tree,
        )
        mocker.patch.object(
            LastCommonAncestorSubtree,
            "read_tree_file",
            side_effect=AssertionError("run should not copy cached trees"),
        )
        mocker.patch(
            "phykit.services.tree.last_common_ancestor_subtree.read_single_column_file_to_list",
            return_value=["a", "b"],
        )
        mocker.patch.object(LastCommonAncestorSubtree, "write_tree_file")

        service.run()

        mocked_read.assert_called_once_with()
        assert tree.common_ancestor_calls == [["a", "b"]]

    def test_run_json_uses_direct_subtree_tip_count(self, mocker, monkeypatch):
        args = Namespace(
            tree="/some/path/to/file.tre",
            list_of_taxa="/some/path/to/taxa.txt",
            output="/tmp/subtree.tre",
            json=True,
        )
        service = LastCommonAncestorSubtree(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        mocker.patch.object(
            LastCommonAncestorSubtree,
            "read_tree_file_unmodified",
            return_value=tree,
        )
        mocker.patch(
            "phykit.services.tree.last_common_ancestor_subtree.read_single_column_file_to_list",
            return_value=["A", "B"],
        )
        mocker.patch.object(LastCommonAncestorSubtree, "write_tree_file")
        mocked_json = mocker.patch(
            "phykit.services.tree.last_common_ancestor_subtree.print_json"
        )

        def fail_count_terminals(*_args, **_kwargs):
            raise AssertionError("JSON count should use direct terminal count")

        monkeypatch.setattr(TreeMixin, "count_terminals", fail_count_terminals)

        service.run()

        assert mocked_json.call_args.args[0]["subtree_tips"] == 2

    def test_find_lca_subtree_avoids_common_ancestor_for_biophylo_tree(self, monkeypatch):
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")

        def fail_common_ancestor(*args, **kwargs):
            raise AssertionError("common_ancestor should not be called")

        monkeypatch.setattr(tree, "common_ancestor", fail_common_ancestor)

        subtree = LastCommonAncestorSubtree._find_lca_subtree(tree, ["A", "B"])

        assert sorted(t.name for t in subtree.get_terminals()) == ["A", "B"]

    def test_find_lca_subtree_uses_direct_traversal(self, monkeypatch):
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")

        def fail_common_ancestor(*args, **kwargs):
            raise AssertionError("common_ancestor should not be called")

        def fail_traversal(*args, **kwargs):
            raise AssertionError("generic traversal should not be called")

        monkeypatch.setattr(tree, "common_ancestor", fail_common_ancestor)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)
        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)

        subtree = LastCommonAncestorSubtree._find_lca_subtree(tree, ["A", "B"])

        assert sorted(t.name for t in subtree.clades) == ["A", "B"]

    def test_find_lca_subtree_single_taxon_uses_direct_terminal_lookup(
        self, monkeypatch
    ):
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")

        def fail_common_ancestor(*args, **kwargs):
            raise AssertionError("single-taxon lookup should not call common_ancestor")

        def fail_traversal(*args, **kwargs):
            raise AssertionError("generic traversal should not be called")

        monkeypatch.setattr(tree, "common_ancestor", fail_common_ancestor)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)
        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)

        subtree = LastCommonAncestorSubtree._find_lca_subtree(tree, ["C", "C"])

        assert subtree.name == "C"

    def test_find_lca_subtree_only_indexes_requested_terminal_names(self, monkeypatch):
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        tree.root.clades[1].clades[0].name = ["unhashable", "non-target"]

        def fail_common_ancestor(*args, **kwargs):
            raise AssertionError("selected-name direct lookup should not fall back")

        monkeypatch.setattr(tree, "common_ancestor", fail_common_ancestor)

        subtree = LastCommonAncestorSubtree._find_lca_subtree(tree, ["A", "B"])

        assert sorted(t.name for t in subtree.clades) == ["A", "B"]

    def test_find_lca_subtree_stops_after_all_requested_taxa_are_found(self):
        left = _Clade(clades=[_Clade(name="A"), _Clade(name="B")])
        root = _Clade(clades=[left, _ExplodingClade()])
        tree = _DirectTree(root)

        subtree = LastCommonAncestorSubtree._find_lca_subtree(tree, ["A", "B"])

        assert subtree is left

    def test_find_lca_subtree_handles_mixed_depth_taxa_with_direct_traversal(
        self, monkeypatch
    ):
        tree = Phylo.read(StringIO("(((A:1,B:1):1,C:1):1,(D:1,E:1):1);"), "newick")

        def fail_common_ancestor(*args, **kwargs):
            raise AssertionError("common_ancestor should not be called")

        monkeypatch.setattr(tree, "common_ancestor", fail_common_ancestor)

        subtree = LastCommonAncestorSubtree._find_lca_subtree(tree, ["A", "C"])

        assert sorted(t.name for t in subtree.get_terminals()) == ["A", "B", "C"]

    def test_find_lca_subtree_handles_multifurcation_with_direct_traversal(
        self, monkeypatch
    ):
        tree = Phylo.read(StringIO("((A:1,B:1,C:1):1,(D:1,E:1):1);"), "newick")

        def fail_common_ancestor(*args, **kwargs):
            raise AssertionError("common_ancestor should not be called")

        def fail_traversal(*args, **kwargs):
            raise AssertionError("generic traversal should not be called")

        monkeypatch.setattr(tree, "common_ancestor", fail_common_ancestor)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)
        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)

        subtree = LastCommonAncestorSubtree._find_lca_subtree(tree, ["A", "C"])

        assert [child.name for child in subtree.clades] == ["A", "B", "C"]
