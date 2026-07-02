from argparse import Namespace
import importlib
import subprocess
import sys
from io import StringIO

import numpy as np
import pytest
import builtins
from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin
from Bio.Phylo.Newick import Clade, Tree

from phykit.services.tree.tip_to_tip_distance import TipToTipDistance


@pytest.fixture
def args():
    return Namespace(tree_zero="/some/path/to/file.tre", tip_1="a", tip_2="b")


class _Tip:
    def __init__(self, name):
        self.name = name


class _Tree:
    def __init__(self, tips):
        self._tips = tips

    def get_terminals(self):
        return [_Tip(tip) for tip in self._tips]


def test_module_import_does_not_import_scipy_clustering(monkeypatch):
    module_name = "phykit.services.tree.tip_to_tip_distance"
    previous = sys.modules.pop(module_name, None)
    original_import = builtins.__import__

    def guarded_import(name, globals=None, locals=None, fromlist=(), level=0):
        if (
            name == "scipy.cluster.hierarchy"
            or name.startswith("scipy.cluster.hierarchy.")
            or name == "scipy.spatial.distance"
            or name.startswith("scipy.spatial.distance.")
        ):
            raise AssertionError(
                "tip_to_tip_distance module import should not import SciPy clustering"
            )
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", guarded_import)
    try:
        importlib.import_module(module_name)
    finally:
        imported = sys.modules.pop(module_name, None)
        if previous is not None:
            sys.modules[module_name] = previous
        parent_name, _, child_name = module_name.rpartition(".")
        parent = sys.modules.get(parent_name)
        if parent is not None:
            if previous is not None:
                setattr(parent, child_name, previous)
            elif getattr(parent, child_name, None) is imported:
                delattr(parent, child_name)


def test_module_import_does_not_import_numpy_or_biopython_tree_modules():
    code = """
import sys
import phykit.services.tree.tip_to_tip_distance as module
assert hasattr(module.np, "__getattr__")
assert hasattr(module.TreeMixin, "distance")
assert "typing" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
assert "numpy" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "Bio.Phylo.BaseTree" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


class TestTipToTipDistance:
    def test_init_sets_expected_attrs(self, args):
        service = TipToTipDistance(args)
        assert service.tree_file_path == args.tree_zero
        assert service.tip_1 == "a"
        assert service.tip_2 == "b"
        assert service.json_output is False
        assert service.all_pairs is False
        assert service.plot is False

    def test_process_args_sets_defaults(self):
        args = Namespace(tree_zero="/some/path/to/file.tre")
        service = TipToTipDistance(args)
        assert service.tip_1 is None
        assert service.tip_2 is None
        assert service.all_pairs is False
        assert service.plot is False
        assert service.plot_output == "tip_to_tip_distance_heatmap.png"

    def test_check_leaves_exits_when_tip_missing(self, mocker, args):
        service = TipToTipDistance(args)
        tree = _Tree(["a", "b"])
        mocker.patch("phykit.services.tree.tip_to_tip_distance.TreeMixin.find_any", side_effect=[object(), None])
        with pytest.raises(SystemExit) as excinfo:
            service.check_leaves(tree, "a", "missing")
        assert excinfo.value.code == 2

    def test_calculate_all_pairwise_distances(self, mocker, args):
        service = TipToTipDistance(args)
        tree = _Tree(["a", "b", "c"])

        def fake_distance(_tree, tip_a, tip_b):
            mapping = {("a", "b"): 1.2345, ("a", "c"): 2.0, ("b", "c"): 3.0}
            return mapping[(tip_a, tip_b)]

        mocker.patch("phykit.services.tree.tip_to_tip_distance.TreeMixin.distance", side_effect=fake_distance)
        rows = service.calculate_all_pairwise_distances(tree)
        assert len(rows) == 3
        assert rows[0] == {"taxon_a": "a", "taxon_b": "b", "tip_to_tip_distance": 1.2345}

    def test_calculate_all_pairwise_distances_fast_matches_biopython(self, args):
        service = TipToTipDistance(args)
        tree = Tree(
            root=Clade(
                clades=[
                    Clade(branch_length=1.0, name="a"),
                    Clade(
                        branch_length=2.0,
                        clades=[
                            Clade(branch_length=3.0, name="b"),
                            Clade(branch_length=4.0, name="c"),
                        ],
                    ),
                ],
            )
        )

        rows = service.calculate_all_pairwise_distances(tree)

        assert rows == [
            {"taxon_a": "a", "taxon_b": "b", "tip_to_tip_distance": round(tree.distance("a", "b"), 4)},
            {"taxon_a": "a", "taxon_b": "c", "tip_to_tip_distance": round(tree.distance("a", "c"), 4)},
            {"taxon_a": "b", "taxon_b": "c", "tip_to_tip_distance": round(tree.distance("b", "c"), 4)},
        ]

    def test_calculate_all_pairwise_distances_uses_direct_tip_names(
        self, monkeypatch, args
    ):
        service = TipToTipDistance(args)
        tree = Tree(
            root=Clade(
                clades=[
                    Clade(branch_length=1.0, name="a"),
                    Clade(branch_length=2.0, name="b"),
                    Clade(branch_length=3.0, name="c"),
                ],
            )
        )

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("all-pairs setup should use direct tip names")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        rows = service.calculate_all_pairwise_distances(tree)

        assert [row["taxon_a"] for row in rows] == ["a", "a", "b"]
        assert [row["taxon_b"] for row in rows] == ["b", "c", "c"]

    def test_calculate_tip_to_tip_distance_fast_matches_biopython(self, args):
        service = TipToTipDistance(args)
        tree = Phylo.read(StringIO("(((A:1,B:2):3,C:4):5,D:6);"), "newick")

        assert service.calculate_tip_to_tip_distance(
            tree, "A", "D"
        ) == pytest.approx(tree.distance("A", "D"))

    def test_calculate_tip_to_tip_distance_fast_uses_child_list(
        self, monkeypatch, args
    ):
        service = TipToTipDistance(args)
        tree = Phylo.read(StringIO("(((A:1,B:2):3,C:4):5,D:6);"), "newick")

        def fail_is_terminal(*_args, **_kwargs):
            raise AssertionError("standard fast path should inspect clades directly")

        monkeypatch.setattr(Clade, "is_terminal", fail_is_terminal)

        assert service.calculate_tip_to_tip_distance(
            tree, "A", "D"
        ) == pytest.approx(15.0)

    def test_calculate_tip_to_tip_distance_binary_path_avoids_reversed(
        self, args
    ):
        class NoReversedList(list):
            def __reversed__(self):
                raise AssertionError("binary fast path should not call reversed")

        def wrap_binary_children(clade):
            children = clade.clades
            if len(children) == 2:
                clade.clades = NoReversedList(children)
            for child in clade.clades:
                wrap_binary_children(child)

        service = TipToTipDistance(args)
        tree = Phylo.read(StringIO("((A:1,B:2):3,(C:4,D:5):6);"), "newick")
        wrap_binary_children(tree.root)

        assert service.calculate_tip_to_tip_distance(
            tree, "A", "D"
        ) == pytest.approx(15.0)

    def test_calculate_tip_to_tip_distance_same_tip_returns_zero_fast(
        self, monkeypatch, args
    ):
        service = TipToTipDistance(args)
        tree = Phylo.read(StringIO("(((A:1,B:2):3,C:4):5,D:6);"), "newick")

        def fail_distance(*_args, **_kwargs):
            raise AssertionError("same-tip standard tree path should return directly")

        monkeypatch.setattr(TreeMixin, "distance", fail_distance)

        assert service.calculate_tip_to_tip_distance(tree, "A", "A") == 0.0

    def test_calculate_tip_to_tip_distance_same_tip_skips_depth_maps(self, args):
        class CladeWithExplodingLength:
            def __init__(self, name=None, clades=None):
                self.name = name
                self.clades = clades or []

            @property
            def branch_length(self):
                raise AssertionError("same-tip lookup should not read branch lengths")

        tree = type(
            "Tree",
            (),
            {
                "root": CladeWithExplodingLength(
                    clades=[
                        CladeWithExplodingLength(
                            clades=[
                                CladeWithExplodingLength(name="A"),
                                CladeWithExplodingLength(name="B"),
                            ]
                        ),
                        CladeWithExplodingLength(name="C"),
                    ]
                )
            },
        )()
        service = TipToTipDistance(args)

        assert service.calculate_tip_to_tip_distance(tree, "C", "C") == 0.0

    def test_build_distance_matrix(self, args):
        service = TipToTipDistance(args)
        rows = [
            {"taxon_a": "a", "taxon_b": "b", "tip_to_tip_distance": 1.0},
            {"taxon_a": "a", "taxon_b": "c", "tip_to_tip_distance": 2.0},
            {"taxon_a": "b", "taxon_b": "c", "tip_to_tip_distance": 3.0},
        ]
        taxa, matrix = service._build_distance_matrix(rows)
        assert taxa == ["a", "b", "c"]
        assert matrix.shape == (3, 3)
        assert np.allclose(matrix, matrix.T)
        assert matrix[0, 1] == 1.0
        assert matrix[1, 2] == 3.0

    def test_build_distance_matrix_uses_sorted_upper_triangle_fast_path(
        self, mocker, monkeypatch, args
    ):
        service = TipToTipDistance(args)
        monkeypatch.setattr(TipToTipDistance, "_MATRIX_FAST_FILL_MIN_ROWS", 0)
        rows = [
            {"taxon_a": "a", "taxon_b": "b", "tip_to_tip_distance": 1.0},
            {"taxon_a": "a", "taxon_b": "c", "tip_to_tip_distance": 2.0},
            {"taxon_a": "b", "taxon_b": "c", "tip_to_tip_distance": 3.0},
        ]
        triu_spy = mocker.spy(
            sys.modules["phykit.services.tree.tip_to_tip_distance"].np,
            "triu_indices",
        )

        taxa, matrix = service._build_distance_matrix(rows)

        assert taxa == ["a", "b", "c"]
        assert triu_spy.call_count == 1
        assert np.allclose(
            matrix,
            np.array(
                [
                    [0.0, 1.0, 2.0],
                    [1.0, 0.0, 3.0],
                    [2.0, 3.0, 0.0],
                ]
            ),
        )

    def test_sorted_upper_triangle_check_skips_unneeded_string_coercion(self):
        class ExplodingStr(str):
            def __str__(self):
                raise AssertionError(
                    "matching row labels should not need string coercion"
                )

        taxa = ["a", "b", "c"]
        rows = [
            {
                "taxon_a": ExplodingStr("a"),
                "taxon_b": ExplodingStr("b"),
                "tip_to_tip_distance": 1.0,
            },
            {
                "taxon_a": ExplodingStr("a"),
                "taxon_b": ExplodingStr("c"),
                "tip_to_tip_distance": 2.0,
            },
            {
                "taxon_a": ExplodingStr("b"),
                "taxon_b": ExplodingStr("c"),
                "tip_to_tip_distance": 3.0,
            },
        ]

        assert TipToTipDistance._rows_are_sorted_upper_triangle(taxa, rows) is True

    def test_sorted_upper_triangle_check_preserves_coercion_fallback(self):
        taxa = ["1", "2"]
        rows = [
            {"taxon_a": 1, "taxon_b": 2, "tip_to_tip_distance": 1.0},
        ]

        assert TipToTipDistance._rows_are_sorted_upper_triangle(taxa, rows) is True

    def test_build_distance_matrix_fallback_handles_arbitrary_row_order(
        self, mocker, args
    ):
        service = TipToTipDistance(args)
        rows = [
            {"taxon_a": "b", "taxon_b": "c", "tip_to_tip_distance": 3.0},
            {"taxon_a": "a", "taxon_b": "c", "tip_to_tip_distance": 2.0},
            {"taxon_a": "a", "taxon_b": "b", "tip_to_tip_distance": 1.0},
        ]
        triu_spy = mocker.spy(
            sys.modules["phykit.services.tree.tip_to_tip_distance"].np,
            "triu_indices",
        )

        taxa, matrix = service._build_distance_matrix(rows)

        assert taxa == ["a", "b", "c"]
        assert triu_spy.call_count == 0
        assert np.allclose(
            matrix,
            np.array(
                [
                    [0.0, 1.0, 2.0],
                    [1.0, 0.0, 3.0],
                    [2.0, 3.0, 0.0],
                ]
            ),
        )

    def test_run_exits_when_plot_without_all_pairs(self, mocker):
        args = Namespace(tree_zero="/some/path/to/file.tre", tip_1="a", tip_2="b", plot=True, all_pairs=False)
        service = TipToTipDistance(args)
        mocker.patch.object(TipToTipDistance, "read_tree_file_unmodified", return_value=_Tree(["a", "b"]))
        with pytest.raises(SystemExit) as excinfo:
            service.run()
        assert excinfo.value.code == 2

    def test_run_pairwise_json(self, mocker):
        args = Namespace(tree_zero="/some/path/to/file.tre", tip_1="a", tip_2="b", json=True)
        service = TipToTipDistance(args)
        mocker.patch.object(TipToTipDistance, "read_tree_file_unmodified", return_value=_Tree(["a", "b"]))
        mocker.patch.object(TipToTipDistance, "check_leaves")
        mocker.patch("phykit.services.tree.tip_to_tip_distance.TreeMixin.distance", return_value=1.23456)
        mocked_json = mocker.patch("phykit.services.tree.tip_to_tip_distance.print_json")
        service.run()
        payload = mocked_json.call_args.args[0]
        assert payload == {"taxon_a": "a", "taxon_b": "b", "tip_to_tip_distance": 1.2346}

    def test_run_pairwise_standard_tree_avoids_distance(self, mocker):
        tree = Phylo.read(StringIO("((A:1,B:2):3,C:4);"), "newick")
        expected = round(tree.distance("A", "C"), 4)
        args = Namespace(tree_zero="/some/path/to/file.tre", tip_1="A", tip_2="C", json=False)
        service = TipToTipDistance(args)
        mocker.patch.object(TipToTipDistance, "read_tree_file_unmodified", return_value=tree)
        mocked_distance = mocker.patch(
            "phykit.services.tree.tip_to_tip_distance.TreeMixin.distance",
            side_effect=AssertionError("standard tree should use fast path"),
        )
        mocked_print = mocker.patch("builtins.print")

        service.run()

        mocked_distance.assert_not_called()
        mocked_print.assert_called_once_with(expected)

    def test_run_all_pairs_json_plot(self, mocker):
        args = Namespace(
            tree_zero="/some/path/to/file.tre",
            all_pairs=True,
            plot=True,
            plot_output="out.png",
            json=True,
            tip_1=None,
            tip_2=None,
        )
        service = TipToTipDistance(args)
        mocker.patch.object(TipToTipDistance, "read_tree_file_unmodified", return_value=_Tree(["a", "b"]))
        mocker.patch.object(
            TipToTipDistance,
            "calculate_all_pairwise_distances",
            return_value=[{"taxon_a": "a", "taxon_b": "b", "tip_to_tip_distance": 1.0}],
        )
        mocked_plot = mocker.patch.object(TipToTipDistance, "_plot_tip_distance_heatmap")
        mocked_json = mocker.patch("phykit.services.tree.tip_to_tip_distance.print_json")
        service.run()
        mocked_plot.assert_called_once()
        payload = mocked_json.call_args.args[0]
        assert payload["all_pairs"] is True
        assert payload["plot_output"] == "out.png"
        assert payload["pairs"] == payload["rows"]

    def test_run_all_pairs_text_output(self, mocker, capsys):
        args = Namespace(
            tree_zero="/some/path/to/file.tre",
            all_pairs=True,
            plot=False,
            json=False,
            tip_1=None,
            tip_2=None,
        )
        service = TipToTipDistance(args)
        mocker.patch.object(TipToTipDistance, "read_tree_file_unmodified", return_value=_Tree(["a", "b"]))
        mocker.patch.object(
            TipToTipDistance,
            "calculate_all_pairwise_distances",
            return_value=[
                {"taxon_a": "a", "taxon_b": "b", "tip_to_tip_distance": 1.0},
                {"taxon_a": "a", "taxon_b": "c", "tip_to_tip_distance": 2.0},
            ],
        )
        service.run()
        captured = capsys.readouterr()
        assert captured.out == "a\tb\t1.0\na\tc\t2.0\n"

    def test_run_all_pairs_text_output_uses_fast_series(self, mocker, capsys):
        args = Namespace(
            tree_zero="/some/path/to/file.tre",
            all_pairs=True,
            plot=False,
            json=False,
            tip_1=None,
            tip_2=None,
        )
        tree = Phylo.read(StringIO("((A:1,B:2):3,C:4);"), "newick")
        service = TipToTipDistance(args)
        mocker.patch.object(TipToTipDistance, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(
            TipToTipDistance,
            "calculate_all_pairwise_distances",
            side_effect=AssertionError("text output should avoid row dictionaries"),
        )

        service.run()

        captured = capsys.readouterr()
        assert captured.out == "A\tB\t3.0\nA\tC\t8.0\nB\tC\t9.0\n"

    def test_run_exits_without_tips_when_not_all_pairs(self, mocker):
        args = Namespace(tree_zero="/some/path/to/file.tre", all_pairs=False, tip_1=None, tip_2=None, json=False)
        service = TipToTipDistance(args)
        mocker.patch.object(TipToTipDistance, "read_tree_file_unmodified", return_value=_Tree(["a", "b"]))
        with pytest.raises(SystemExit) as excinfo:
            service.run()
        assert excinfo.value.code == 2

    def test_run_pairwise_text_output(self, mocker, capsys):
        args = Namespace(tree_zero="/some/path/to/file.tre", tip_1="a", tip_2="b", json=False, all_pairs=False)
        service = TipToTipDistance(args)
        mocker.patch.object(TipToTipDistance, "read_tree_file_unmodified", return_value=_Tree(["a", "b"]))
        mocker.patch.object(TipToTipDistance, "check_leaves")
        mocker.patch("phykit.services.tree.tip_to_tip_distance.TreeMixin.distance", return_value=2.34567)
        service.run()
        out, _ = capsys.readouterr()
        assert out.strip() == "2.3457"

    def test_run_uses_unmodified_tree_read(self, mocker):
        tree = _Tree(["a", "b"])
        args = Namespace(tree_zero="/some/path/to/file.tre", tip_1="a", tip_2="b", json=False, all_pairs=False, plot=False)
        service = TipToTipDistance(args)
        read_tree = mocker.patch.object(service, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(
            service,
            "read_tree_file",
            side_effect=AssertionError("copying tree reader should not be used"),
        )
        calc = mocker.patch.object(service, "calculate_tip_to_tip_distance", return_value=2.34567)

        service.run()

        read_tree.assert_called_once_with()
        calc.assert_called_once_with(tree, "a", "b")

    def test_plot_tip_distance_heatmap_creates_file(self, tmp_path):
        pytest.importorskip("matplotlib")
        out = tmp_path / "tip_heatmap.png"
        service = TipToTipDistance(
            Namespace(tree_zero="/some/path/to/file.tre", all_pairs=True, plot=True, plot_output=str(out))
        )
        service._plot_tip_distance_heatmap(
            [
                {"taxon_a": "a", "taxon_b": "b", "tip_to_tip_distance": 1.0},
                {"taxon_a": "a", "taxon_b": "c", "tip_to_tip_distance": 2.0},
                {"taxon_a": "b", "taxon_b": "c", "tip_to_tip_distance": 3.0},
            ]
        )
        assert out.exists()

    def test_plot_tip_distance_heatmap_empty_rows(self):
        pytest.importorskip("matplotlib")
        service = TipToTipDistance(Namespace(tree_zero="/some/path/to/file.tre", all_pairs=True, plot=True))
        service._plot_tip_distance_heatmap([])

    def test_plot_tip_distance_heatmap_importerror(self, monkeypatch, capsys):
        original_import = builtins.__import__

        def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name.startswith("matplotlib"):
                raise ImportError("no matplotlib")
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", fake_import)
        service = TipToTipDistance(Namespace(tree_zero="/some/path/to/file.tre", all_pairs=True, plot=True))
        with pytest.raises(SystemExit) as exc:
            service._plot_tip_distance_heatmap([{"taxon_a": "a", "taxon_b": "b", "tip_to_tip_distance": 1.0}])
        assert exc.value.code == 2
        out, _ = capsys.readouterr()
        assert "matplotlib is required for --plot in tip_to_tip_distance" in out

    def test_check_leaves_first_missing(self, mocker, args):
        service = TipToTipDistance(args)
        tree = _Tree(["a", "b"])
        mocker.patch("phykit.services.tree.tip_to_tip_distance.TreeMixin.find_any", return_value=None)
        with pytest.raises(SystemExit) as excinfo:
            service.check_leaves(tree, "missing", "b")
        assert excinfo.value.code == 2
