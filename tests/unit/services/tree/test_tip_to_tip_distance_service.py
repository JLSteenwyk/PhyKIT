from argparse import Namespace

import numpy as np
import pytest
import builtins

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

    def test_run_exits_when_plot_without_all_pairs(self, mocker):
        args = Namespace(tree_zero="/some/path/to/file.tre", tip_1="a", tip_2="b", plot=True, all_pairs=False)
        service = TipToTipDistance(args)
        mocker.patch.object(TipToTipDistance, "read_tree_file", return_value=_Tree(["a", "b"]))
        with pytest.raises(SystemExit) as excinfo:
            service.run()
        assert excinfo.value.code == 2

    def test_run_pairwise_json(self, mocker):
        args = Namespace(tree_zero="/some/path/to/file.tre", tip_1="a", tip_2="b", json=True)
        service = TipToTipDistance(args)
        mocker.patch.object(TipToTipDistance, "read_tree_file", return_value=_Tree(["a", "b"]))
        mocker.patch.object(TipToTipDistance, "check_leaves")
        mocker.patch("phykit.services.tree.tip_to_tip_distance.TreeMixin.distance", return_value=1.23456)
        mocked_json = mocker.patch("phykit.services.tree.tip_to_tip_distance.print_json")
        service.run()
        payload = mocked_json.call_args.args[0]
        assert payload == {"taxon_a": "a", "taxon_b": "b", "tip_to_tip_distance": 1.2346}

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
        mocker.patch.object(TipToTipDistance, "read_tree_file", return_value=_Tree(["a", "b"]))
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

    def test_run_all_pairs_text_output(self, mocker):
        args = Namespace(
            tree_zero="/some/path/to/file.tre",
            all_pairs=True,
            plot=False,
            json=False,
            tip_1=None,
            tip_2=None,
        )
        service = TipToTipDistance(args)
        mocker.patch.object(TipToTipDistance, "read_tree_file", return_value=_Tree(["a", "b"]))
        mocker.patch.object(
            TipToTipDistance,
            "calculate_all_pairwise_distances",
            return_value=[{"taxon_a": "a", "taxon_b": "b", "tip_to_tip_distance": 1.0}],
        )
        mocked_print = mocker.patch("builtins.print")
        service.run()
        mocked_print.assert_called_once_with("a\tb\t1.0")

    def test_run_exits_without_tips_when_not_all_pairs(self, mocker):
        args = Namespace(tree_zero="/some/path/to/file.tre", all_pairs=False, tip_1=None, tip_2=None, json=False)
        service = TipToTipDistance(args)
        mocker.patch.object(TipToTipDistance, "read_tree_file", return_value=_Tree(["a", "b"]))
        with pytest.raises(SystemExit) as excinfo:
            service.run()
        assert excinfo.value.code == 2

    def test_run_pairwise_text_output(self, mocker, capsys):
        args = Namespace(tree_zero="/some/path/to/file.tre", tip_1="a", tip_2="b", json=False, all_pairs=False)
        service = TipToTipDistance(args)
        mocker.patch.object(TipToTipDistance, "read_tree_file", return_value=_Tree(["a", "b"]))
        mocker.patch.object(TipToTipDistance, "check_leaves")
        mocker.patch("phykit.services.tree.tip_to_tip_distance.TreeMixin.distance", return_value=2.34567)
        service.run()
        out, _ = capsys.readouterr()
        assert out.strip() == "2.3457"

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
