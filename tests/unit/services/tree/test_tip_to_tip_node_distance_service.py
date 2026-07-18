from argparse import Namespace
from io import StringIO
import subprocess
import sys

import pytest
from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin

from phykit.services.tree.tip_to_tip_node_distance import TipToTipNodeDistance
import phykit.services.tree.tip_to_tip_node_distance as node_distance_module


def test_module_import_does_not_import_biophylo_or_numpy():
    code = """
import sys
import phykit.services.tree.tip_to_tip_node_distance as module

assert hasattr(module.TreeMixin, "find_any")
assert hasattr(module.TreeMixin, "trace")
assert "json" not in sys.modules
assert "typing" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "numpy" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_lazy_tree_mixin_caches_resolved_find_any_and_trace(monkeypatch):
    lazy_tree_mixin = node_distance_module._LazyTreeMixin()

    monkeypatch.setattr(TreeMixin, "find_any", lambda *args, **kwargs: "found")
    monkeypatch.setattr(TreeMixin, "trace", lambda *args, **kwargs: [1, 2, 3])

    assert lazy_tree_mixin.find_any("tree", "A") == "found"
    assert lazy_tree_mixin.trace("tree", "A", "B") == [1, 2, 3]
    cached_find_any = lazy_tree_mixin.__dict__["find_any"]
    cached_trace = lazy_tree_mixin.__dict__["trace"]

    monkeypatch.setattr(
        TreeMixin,
        "find_any",
        lambda *args, **kwargs: pytest.fail("cached find_any should be reused"),
    )
    monkeypatch.setattr(
        TreeMixin,
        "trace",
        lambda *args, **kwargs: pytest.fail("cached trace should be reused"),
    )

    assert lazy_tree_mixin.find_any("tree", "A") == "found"
    assert lazy_tree_mixin.trace("tree", "A", "B") == [1, 2, 3]
    assert lazy_tree_mixin.__dict__["find_any"] is cached_find_any
    assert lazy_tree_mixin.__dict__["trace"] is cached_trace


@pytest.fixture
def args():
    return Namespace(tree_zero="/some/path/to/file.tre", tip_1="a", tip_2="b")


class _Tree:
    pass


class _ChildListClade:
    def __init__(self, name=None, clades=None):
        self.name = name
        self.clades = clades or []

    def is_terminal(self):
        raise AssertionError("fast path should inspect clades directly")


class TestTipToTipNodeDistance:
    def test_init_sets_expected_attrs(self, args):
        service = TipToTipNodeDistance(args)
        assert service.tree_file_path == args.tree_zero
        assert service.tip_1 == "a"
        assert service.tip_2 == "b"
        assert service.json_output is False

    def test_check_leaves_exits_when_first_missing(self, mocker, args):
        service = TipToTipNodeDistance(args)
        mocker.patch(
            "phykit.services.tree.tip_to_tip_node_distance.TreeMixin.find_any",
            side_effect=[None, object()],
        )
        with pytest.raises(SystemExit) as excinfo:
            service.check_leaves(_Tree(), "missing", "b")
        assert excinfo.value.code == 2

    def test_check_leaves_exits_when_second_missing(self, mocker, args):
        service = TipToTipNodeDistance(args)
        mocker.patch(
            "phykit.services.tree.tip_to_tip_node_distance.TreeMixin.find_any",
            side_effect=[object(), None],
        )

        with pytest.raises(SystemExit) as excinfo:
            service.check_leaves(_Tree(), "a", "missing")

        assert excinfo.value.code == 2

    def test_run_text_output(self, mocker):
        args = Namespace(tree_zero="/some/path/to/file.tre", tip_1="a", tip_2="b", json=False)
        service = TipToTipNodeDistance(args)
        mocker.patch.object(TipToTipNodeDistance, "read_tree_file_unmodified", return_value=_Tree())
        mocker.patch.object(TipToTipNodeDistance, "check_leaves")
        mocker.patch("phykit.services.tree.tip_to_tip_node_distance.TreeMixin.trace", return_value=[1, 2, 3, 4])
        mocked_print = mocker.patch("builtins.print")
        service.run()
        mocked_print.assert_called_once_with(4)

    def test_run_json_output(self, mocker):
        args = Namespace(tree_zero="/some/path/to/file.tre", tip_1="a", tip_2="b", json=True)
        service = TipToTipNodeDistance(args)
        mocker.patch.object(TipToTipNodeDistance, "read_tree_file_unmodified", return_value=_Tree())
        mocker.patch.object(TipToTipNodeDistance, "check_leaves")
        mocker.patch("phykit.services.tree.tip_to_tip_node_distance.TreeMixin.trace", return_value=[1, 2, 3])
        mocked_json = mocker.patch("phykit.services.tree.tip_to_tip_node_distance.print_json")
        service.run()
        payload = mocked_json.call_args.args[0]
        assert payload == {
            "taxon_a": "a",
            "taxon_b": "b",
            "tip_to_tip_node_distance": 3,
        }

    def test_fast_node_distance_matches_biopython_trace(self, args):
        tree = Phylo.read(StringIO("(((A:1,B:1):1,C:2):1,D:3);"), "newick")
        service = TipToTipNodeDistance(args)

        assert service.calculate_tip_to_tip_node_distance(
            tree, "A", "D"
        ) == len(TreeMixin.trace(tree, "A", "D"))

    def test_fast_node_distance_uses_child_lists_for_terminal_checks(self, args):
        leaf_a = _ChildListClade("A")
        leaf_b = _ChildListClade("B")
        leaf_c = _ChildListClade("C")
        internal = _ChildListClade(clades=[leaf_a, leaf_b])
        tree = _Tree()
        tree.root = _ChildListClade(clades=[internal, leaf_c])
        service = TipToTipNodeDistance(args)

        assert service.calculate_tip_to_tip_node_distance(tree, "A", "C") == 3

    def test_fast_node_distance_supports_multifurcations(self, args):
        tree = Phylo.read(StringIO("(A:1,B:1,C:1,D:1);"), "newick")
        service = TipToTipNodeDistance(args)

        observed = service.calculate_tip_to_tip_node_distance(tree, "A", "D")

        assert observed == len(TreeMixin.trace(tree, "A", "D"))

    def test_fast_node_distance_reports_missing_first_tip(self, args, capsys):
        tree = Phylo.read(StringIO("((A:1,B:1):1,C:1);"), "newick")
        service = TipToTipNodeDistance(args)

        with pytest.raises(SystemExit) as excinfo:
            service.calculate_tip_to_tip_node_distance(tree, "missing", "B")

        assert excinfo.value.code == 2
        assert "missing not on tree" in capsys.readouterr().out

    def test_fast_node_distance_does_not_build_root_paths(self, args, mocker):
        tree = Phylo.read(StringIO("(((A:1,B:1):1,C:2):1,D:3);"), "newick")
        service = TipToTipNodeDistance(args)
        mocker.patch.object(
            TipToTipNodeDistance,
            "_path_to_root",
            side_effect=AssertionError("fast path should use ancestor distances"),
        )

        assert service.calculate_tip_to_tip_node_distance(tree, "A", "D") == 4

    def test_fast_node_distance_same_tip_returns_zero_without_trace(
        self,
        args,
        mocker,
    ):
        tree = Phylo.read(StringIO("(((A:1,B:1):1,C:2):1,D:3);"), "newick")
        service = TipToTipNodeDistance(args)
        mocker.patch(
            "phykit.services.tree.tip_to_tip_node_distance.TreeMixin.trace",
            side_effect=AssertionError("same-tip standard tree path should return directly"),
        )

        assert service.calculate_tip_to_tip_node_distance(tree, "A", "A") == 0

    def test_fast_node_distance_same_tip_skips_parent_map(self, args, mocker):
        leaf_a = _ChildListClade("A")
        leaf_b = _ChildListClade("B")
        leaf_c = _ChildListClade("C")
        internal = _ChildListClade(clades=[leaf_a, leaf_b])
        tree = _Tree()
        tree.root = _ChildListClade(clades=[internal, leaf_c])
        service = TipToTipNodeDistance(args)
        mocker.patch.object(
            TipToTipNodeDistance,
            "_path_to_root",
            side_effect=AssertionError("same-tip path should not build root paths"),
        )

        assert service.calculate_tip_to_tip_node_distance(tree, "C", "C") == 0

    def test_fast_same_tip_supports_multifurcations(self, args):
        tree = Phylo.read(StringIO("(A:1,B:1,C:1,D:1);"), "newick")
        service = TipToTipNodeDistance(args)

        assert service.calculate_tip_to_tip_node_distance(tree, "D", "D") == 0

    def test_fast_same_tip_returns_none_for_incomplete_clade_interface(self):
        root = _ChildListClade(clades=[object()])

        observed = TipToTipNodeDistance._calculate_same_tip_node_distance_fast(
            root, "A"
        )

        assert observed is None

    def test_fast_same_tip_reports_missing_taxon(self, args, capsys):
        tree = Phylo.read(StringIO("((A:1,B:1):1,C:1);"), "newick")
        service = TipToTipNodeDistance(args)

        with pytest.raises(SystemExit) as excinfo:
            service.calculate_tip_to_tip_node_distance(
                tree, "missing", "missing"
            )

        assert excinfo.value.code == 2
        assert "missing not on tree" in capsys.readouterr().out

    def test_path_to_root_returns_nodes_in_root_first_order(self, args):
        tree = Phylo.read(StringIO("((A:1,B:1):1,C:1);"), "newick")
        service = TipToTipNodeDistance(args)
        leaf = next(tree.find_clades(name="A"))
        parent = tree.root.clades[0]
        parent_map = {leaf: parent, parent: tree.root}

        path = service._path_to_root(leaf, tree.root, parent_map)

        assert path == [tree.root, parent, leaf]

    def test_run_standard_tree_avoids_trace(self, mocker):
        tree = Phylo.read(StringIO("((A:1,B:1):1,C:2);"), "newick")
        expected = len(TreeMixin.trace(tree, "A", "C"))
        args = Namespace(tree_zero="/some/path/to/file.tre", tip_1="A", tip_2="C", json=False)
        service = TipToTipNodeDistance(args)
        mocker.patch.object(TipToTipNodeDistance, "read_tree_file_unmodified", return_value=tree)
        mocked_trace = mocker.patch(
            "phykit.services.tree.tip_to_tip_node_distance.TreeMixin.trace",
            side_effect=AssertionError("standard tree should use fast path"),
        )
        mocked_print = mocker.patch("builtins.print")

        service.run()

        mocked_trace.assert_not_called()
        mocked_print.assert_called_once_with(expected)

    def test_run_uses_unmodified_tree_read(self, mocker):
        tree = _Tree()
        args = Namespace(tree_zero="/some/path/to/file.tre", tip_1="a", tip_2="b", json=False)
        service = TipToTipNodeDistance(args)
        read_tree = mocker.patch.object(
            service, "read_tree_file_unmodified", return_value=tree
        )
        mocker.patch.object(
            service,
            "read_tree_file",
            side_effect=AssertionError("copying tree reader should not be used"),
        )
        calc = mocker.patch.object(
            service, "calculate_tip_to_tip_node_distance", return_value=3
        )
        mocker.patch("builtins.print")

        service.run()

        read_tree.assert_called_once_with()
        calc.assert_called_once_with(tree, "a", "b")
