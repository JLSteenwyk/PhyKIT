import pytest
from argparse import Namespace
from io import StringIO
import subprocess
import sys

from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin
from phykit.services.tree.print_tree import PrintTree
import phykit.services.tree.print_tree as print_tree_module


@pytest.fixture
def args():
    kwargs = dict(tree="/some/path/to/file.tre", remove=None)
    return Namespace(**kwargs)


def test_module_import_does_not_import_biophylo_or_numpy():
    code = """
import sys
import phykit.services.tree.print_tree as module

assert hasattr(module.Phylo, "draw_ascii")
assert "typing" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "numpy" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_lazy_phylo_caches_resolved_draw_ascii(monkeypatch):
    calls = []

    def cached_draw_ascii(*args, **kwargs):
        calls.append((args, kwargs))
        return "cached"

    def uncached_draw_ascii(*_args, **_kwargs):
        return "uncached"

    lazy_phylo = print_tree_module._LazyPhylo()
    monkeypatch.setattr(Phylo, "draw_ascii", cached_draw_ascii)

    assert lazy_phylo.draw_ascii("tree") == "cached"

    monkeypatch.setattr(Phylo, "draw_ascii", uncached_draw_ascii)

    assert lazy_phylo.draw_ascii("tree2") == "cached"
    assert lazy_phylo.__dict__["draw_ascii"] is cached_draw_ascii
    assert calls == [(("tree",), {}), (("tree2",), {})]


class TestPrintTree(object):
    def test_init_sets_tree_file_path(self, args):
        t = PrintTree(args)
        assert t.tree_file_path == args.tree
        assert t.output_file_path is None

    def test_read_file_reads_tree_file_path(self, mocker, args):
        # Mock the cached tree read method instead of Phylo.read
        mock_cached_read = mocker.patch("phykit.services.tree.base.Tree._cached_tree_read")
        mock_get_hash = mocker.patch("phykit.services.tree.base.Tree._get_file_hash", return_value="test_hash")

        t = PrintTree(args)
        t.read_tree_file()

        # Verify the cached read was called with the correct parameters
        mock_get_hash.assert_called_with(args.tree)
        mock_cached_read.assert_called_with(args.tree, "newick", "test_hash")

    def test_process_args_defaults_json_false(self):
        parsed = PrintTree(Namespace(tree="t.tre", remove=False)).process_args(
            Namespace(tree="t.tre", remove=False)
        )
        assert parsed["json_output"] is False

    def test_run_json_remove_branch_lengths(self, tree_simple, monkeypatch):
        captured = {}
        before_lengths = [node.branch_length for node in tree_simple.get_terminals()]
        t = PrintTree(Namespace(tree="x.tre", remove=True, json=True))
        monkeypatch.setattr(t, "read_tree_file", lambda: t._fast_copy(tree_simple))
        monkeypatch.setattr(print_tree_module, "print_json", lambda payload: captured.setdefault("payload", payload))

        t.run()

        assert captured["payload"]["remove_branch_lengths"] is True
        assert "0.00000" in captured["payload"]["tree_newick"]
        assert [node.branch_length for node in tree_simple.get_terminals()] == before_lengths

    def test_remove_branch_lengths_uses_direct_traversal(self, args, monkeypatch):
        tree = Phylo.read(StringIO("((A:2,B:3):5,C:7);"), "newick")
        t = PrintTree(args)

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("generic tree traversal should not be used")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        t.remove_branch_lengths(tree)

        assert tree.root.branch_length is None
        assert tree.root.clades[0].branch_length is None
        assert tree.root.clades[0].clades[0].branch_length is None
        assert tree.root.clades[0].clades[1].branch_length is None
        assert tree.root.clades[1].branch_length is None

    def test_remove_standard_tree_branch_lengths_clears_during_direct_traversal(self):
        tree = Phylo.read(StringIO("((A:2,B:3):5,C:7);"), "newick")

        assert PrintTree._remove_standard_tree_branch_lengths(tree) is True

        assert all(node.branch_length is None for node in tree.find_clades())

    def test_remove_standard_tree_branch_lengths_handles_mixed_child_counts(self):
        class Clade:
            def __init__(self, branch_length=None, clades=None):
                self.branch_length = branch_length
                self.clades = clades or []

        tree = type(
            "Tree",
            (),
            {
                "root": Clade(
                    0.0,
                    [
                        Clade(1.0),
                        Clade(2.0, [Clade(3.0)]),
                        Clade(4.0, [Clade(5.0), Clade(None), Clade(6.0)]),
                    ],
                )
            },
        )()

        assert PrintTree._remove_standard_tree_branch_lengths(tree) is True

        stack = [tree.root]
        while stack:
            clade = stack.pop()
            assert clade.branch_length is None
            stack.extend(clade.clades)

    def test_run_draw_ascii_and_broken_pipe(self, tree_simple, monkeypatch):
        t = PrintTree(Namespace(tree="x.tre", remove=False, json=False))
        monkeypatch.setattr(t, "read_tree_file_unmodified", lambda: tree_simple)
        monkeypatch.setattr(print_tree_module.Phylo, "draw_ascii", lambda _tree: (_ for _ in ()).throw(BrokenPipeError()))

        # no exception should be raised
        t.run()

    def test_run_no_remove_uses_unmodified_tree_read(self, mocker):
        tree = object()
        t = PrintTree(Namespace(tree="x.tre", remove=False, json=False))
        read_tree = mocker.patch.object(t, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(
            t,
            "read_tree_file",
            side_effect=AssertionError("copying tree reader should not be used"),
        )
        draw_ascii = mocker.patch.object(print_tree_module.Phylo, "draw_ascii")

        t.run()

        read_tree.assert_called_once_with()
        draw_ascii.assert_called_once_with(tree)

    def test_run_remove_uses_copying_tree_read(self, tree_simple, mocker):
        t = PrintTree(Namespace(tree="x.tre", remove=True, json=False))
        read_tree = mocker.patch.object(t, "read_tree_file", return_value=tree_simple)
        mocker.patch.object(
            t,
            "read_tree_file_unmodified",
            side_effect=AssertionError("unmodified reader should not be used"),
        )
        mocker.patch.object(print_tree_module.Phylo, "draw_ascii")

        t.run()

        read_tree.assert_called_once_with()
