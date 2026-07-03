from argparse import Namespace
from io import StringIO
import json
import subprocess
import sys

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree, TreeMixin
import pytest

from phykit.helpers.plot_config import build_parent_map
import phykit.services.tree.chronogram as module
from phykit.services.tree.chronogram import Chronogram


TREE = "tests/sample_files/ultrametric_tree.tre"


def test_module_import_does_not_import_numpy_matplotlib_or_biophylo():
    code = """
import sys
import phykit.services.tree.chronogram as module

assert hasattr(module.np, "__getattr__")
assert callable(module.print_json)
assert callable(module.build_parent_map)
assert callable(module.compute_node_positions)
assert callable(module.compute_circular_coords)
assert "typing" not in sys.modules
assert "json" not in sys.modules
assert "numpy" not in sys.modules
assert "matplotlib" not in sys.modules
assert "matplotlib.pyplot" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
assert "phykit.helpers.geological_timescale" not in sys.modules
assert "phykit.helpers.circular_layout" not in sys.modules
assert "phykit.helpers.color_annotations" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_lazy_numpy_caches_resolved_attributes():
    lazy_np = module._LazyNumpy()

    linspace_attr = lazy_np.linspace

    assert lazy_np.__dict__["linspace"] is linspace_attr
    assert lazy_np.linspace is linspace_attr
    assert lazy_np._module is not None


def _make_args(**overrides):
    defaults = dict(
        tree=TREE, root_age=70.0, plot_output=None, timescale="auto",
        node_ages=False, json=False,
        fig_width=None, fig_height=None, dpi=150, no_title=False,
        title=None, legend_position=None, ylabel_fontsize=None,
        xlabel_fontsize=None, title_fontsize=None, axis_fontsize=None,
        colors=None, ladderize=False, cladogram=False, circular=False,
        color_file=None,
    )
    defaults.update(overrides)
    return Namespace(**defaults)


class TestChronogram:
    def test_run_reuses_unmodified_tree(self, mocker, tmp_path):
        out = str(tmp_path / "chrono.png")
        svc = Chronogram(_make_args(plot_output=out))
        tree = object()
        root_to_tip = {1: 0.0, 2: 1.0}

        mocker.patch.object(
            svc,
            "read_tree_file",
            side_effect=AssertionError("run should not copy the cached tree"),
        )
        read_unmodified = mocker.patch.object(
            svc, "read_tree_file_unmodified", return_value=tree
        )
        validate = mocker.patch.object(svc, "validate_tree")
        parent_map = mocker.patch(
            "phykit.services.tree.chronogram.build_parent_map",
            return_value={},
        )
        compute = mocker.patch.object(
            svc, "_compute_root_to_tip", return_value=root_to_tip
        )
        hpd = mocker.patch.object(svc, "_parse_hpd_intervals", return_value={})
        plot = mocker.patch.object(svc, "_plot_rectangular")

        svc.run()

        read_unmodified.assert_called_once_with()
        validate.assert_called_once_with(
            tree, min_tips=3, require_branch_lengths=True, context="chronogram"
        )
        parent_map.assert_called_once_with(tree)
        compute.assert_called_once_with(tree)
        hpd.assert_called_once_with(tree, 70.0, 70.0)
        plot.assert_called_once_with(tree, {}, 70.0, root_to_tip, {})

    def test_ladderized_run_copies_before_mutating_tree(self, mocker, tmp_path):
        out = str(tmp_path / "chrono_ladder.png")
        svc = Chronogram(_make_args(plot_output=out, ladderize=True))
        cached_tree = mocker.Mock()
        plot_tree = mocker.Mock()
        root_to_tip = {1: 0.0, 2: 1.0}

        mocker.patch.object(svc, "read_tree_file_unmodified", return_value=cached_tree)
        copy_tree = mocker.patch.object(svc, "_fast_copy", return_value=plot_tree)
        mocker.patch.object(svc, "validate_tree")
        mocker.patch(
            "phykit.services.tree.chronogram.build_parent_map",
            return_value={},
        )
        mocker.patch.object(svc, "_compute_root_to_tip", return_value=root_to_tip)
        mocker.patch.object(svc, "_parse_hpd_intervals", return_value={})
        plot = mocker.patch.object(svc, "_plot_rectangular")

        svc.run()

        copy_tree.assert_called_once_with(cached_tree)
        cached_tree.ladderize.assert_not_called()
        plot_tree.ladderize.assert_called_once_with()
        plot.assert_called_once()
        assert plot.call_args.args[0] is plot_tree

    def test_parse_hpd_intervals_uses_direct_tree_traversal(self, monkeypatch):
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        tree.root.clades[0].comment = "[&height_95%_HPD={1.25,2.5}]"
        tree.root.clades[1].comment = "[&95%HPD={3.0,1.5}]"

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("standard tree should use direct preorder traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        intervals = Chronogram._parse_hpd_intervals(tree, 70.0, 1.0)

        assert intervals[id(tree.root.clades[0])] == (1.25, 2.5)
        assert intervals[id(tree.root.clades[1])] == (1.5, 3.0)

    def test_parse_hpd_intervals_ignores_missing_comment_attribute(self):
        tree = Tree(root=Clade(clades=[Clade(name="A"), Clade(name="B")]))

        intervals = Chronogram._parse_hpd_intervals(tree, 70.0, 1.0)

        assert intervals == {}

    @pytest.mark.parametrize("circular", [False, True])
    def test_plot_uses_direct_tree_traversal(self, monkeypatch, tmp_path, circular):
        out = str(tmp_path / f"chrono_direct_{circular}.png")
        svc = Chronogram(_make_args(
            plot_output=out, circular=circular, node_ages=True,
        ))
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        parent_map = build_parent_map(tree)
        root_to_tip = svc._compute_root_to_tip(tree)
        hpd_intervals = {
            id(clade): (10.0, 20.0)
            for clade in svc._iter_preorder(tree.root)
            if clade.clades and clade != tree.root
        }

        def fail_traversal(*args, **kwargs):
            raise AssertionError("plot setup should use direct traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)

        if circular:
            svc._plot_circular(tree, parent_map, 1.0, root_to_tip, hpd_intervals)
        else:
            svc._plot_rectangular(tree, parent_map, 1.0, root_to_tip, hpd_intervals)

        assert (tmp_path / f"chrono_direct_{circular}.png").exists()

    @pytest.mark.parametrize("circular", [False, True])
    def test_plot_reuses_clade_lists_for_layout_helpers(
        self, monkeypatch, tmp_path, circular
    ):
        out = str(tmp_path / f"chrono_layout_lists_{circular}.png")
        svc = Chronogram(_make_args(plot_output=out, circular=circular))
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        parent_map = build_parent_map(tree)
        root_to_tip = svc._compute_root_to_tip(tree)
        layout_calls = {}
        original_node_positions = module.compute_node_positions
        original_circular_coords = module.compute_circular_coords

        def spy_node_positions(*args, **kwargs):
            layout_calls["node_preorder"] = kwargs.get("preorder_clades")
            return original_node_positions(*args, **kwargs)

        def spy_circular_coords(*args, **kwargs):
            layout_calls["circular_preorder"] = kwargs.get("preorder_clades")
            layout_calls["circular_tips"] = kwargs.get("terminal_clades")
            return original_circular_coords(*args, **kwargs)

        monkeypatch.setattr(module, "compute_node_positions", spy_node_positions)
        monkeypatch.setattr(module, "compute_circular_coords", spy_circular_coords)

        if circular:
            svc._plot_circular(tree, parent_map, 1.0, root_to_tip)
        else:
            svc._plot_rectangular(tree, parent_map, 1.0, root_to_tip)

        expected_preorder = list(svc._iter_preorder(tree.root))
        if circular:
            expected_tips = [
                clade for clade in expected_preorder if not clade.clades
            ]
            assert "node_preorder" not in layout_calls
            assert len(layout_calls["circular_preorder"]) == len(expected_preorder)
            assert len(layout_calls["circular_tips"]) == len(expected_tips)
            assert all(
                actual is expected
                for actual, expected in zip(
                    layout_calls["circular_preorder"], expected_preorder
                )
            )
            assert all(
                actual is expected
                for actual, expected in zip(
                    layout_calls["circular_tips"], expected_tips
                )
            )
        else:
            assert "circular_preorder" not in layout_calls
            assert len(layout_calls["node_preorder"]) == len(expected_preorder)
            assert all(
                actual is expected
                for actual, expected in zip(
                    layout_calls["node_preorder"], expected_preorder
                )
            )
        assert (tmp_path / f"chrono_layout_lists_{circular}.png").exists()

    def test_rectangular_plot_batches_base_branches(self, monkeypatch, tmp_path):
        pytest.importorskip("matplotlib")
        import matplotlib.axes

        out = str(tmp_path / "chrono_batched.png")
        svc = Chronogram(_make_args(plot_output=out))
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        parent_map = build_parent_map(tree)
        root_to_tip = svc._compute_root_to_tip(tree)

        def fail_plot(*args, **kwargs):
            raise AssertionError("rectangular chronogram branches should be batched")

        monkeypatch.setattr(matplotlib.axes.Axes, "plot", fail_plot)

        svc._plot_rectangular(tree, parent_map, 1.0, root_to_tip)

        assert (tmp_path / "chrono_batched.png").exists()

    def test_rectangular_plot_batches_hpd_intervals(self, monkeypatch, tmp_path):
        pytest.importorskip("matplotlib")
        import matplotlib.axes
        from matplotlib.collections import PatchCollection

        out = str(tmp_path / "chrono_rectangular_hpd_batched.png")
        svc = Chronogram(_make_args(plot_output=out))
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        parent_map = build_parent_map(tree)
        root_to_tip = svc._compute_root_to_tip(tree)
        hpd_intervals = {
            id(clade): (10.0, 20.0)
            for clade in svc._iter_preorder(tree.root)
            if clade.clades and clade != tree.root
        }
        original_add_collection = matplotlib.axes.Axes.add_collection
        patch_collections = []

        def fail_barh(*args, **kwargs):
            raise AssertionError("rectangular HPD intervals should be batched")

        def count_collection(self, collection, *args, **kwargs):
            if isinstance(collection, PatchCollection):
                patch_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "barh", fail_barh)
        monkeypatch.setattr(
            matplotlib.axes.Axes, "add_collection", count_collection
        )

        svc._plot_rectangular(tree, parent_map, 1.0, root_to_tip, hpd_intervals)

        assert any(
            len(collection.get_paths()) == len(hpd_intervals)
            for collection in patch_collections
        )
        assert (tmp_path / "chrono_rectangular_hpd_batched.png").exists()

    def test_rectangular_plot_batches_color_clade_overlay(
        self, monkeypatch, tmp_path
    ):
        pytest.importorskip("matplotlib")
        import matplotlib.axes
        from matplotlib.collections import LineCollection

        color_file = tmp_path / "colors.tsv"
        color_file.write_text("A,B\tclade\t#ff0000\tAB\n")
        out = str(tmp_path / "chrono_color_batched.png")
        svc = Chronogram(_make_args(plot_output=out, color_file=str(color_file)))
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        parent_map = build_parent_map(tree)
        root_to_tip = svc._compute_root_to_tip(tree)
        original_plot = matplotlib.axes.Axes.plot
        original_add_collection = matplotlib.axes.Axes.add_collection
        line_collections = []

        def fail_overlay_plot(self, *args, **kwargs):
            if kwargs.get("color") == "#ff0000" and kwargs.get("lw") == 1.5:
                raise AssertionError("clade overlay branches should be batched")
            return original_plot(self, *args, **kwargs)

        def count_collection(self, collection, *args, **kwargs):
            if isinstance(collection, LineCollection):
                line_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "plot", fail_overlay_plot)
        monkeypatch.setattr(
            matplotlib.axes.Axes, "add_collection", count_collection
        )

        svc._plot_rectangular(tree, parent_map, 1.0, root_to_tip)

        assert any(
            len(collection.get_segments()) >= 1
            and collection.get_linewidths()[0] == 1.5
            for collection in line_collections
        )
        assert (tmp_path / "chrono_color_batched.png").exists()

    def test_circular_plot_batches_base_branches(self, monkeypatch, tmp_path):
        pytest.importorskip("matplotlib")
        import matplotlib.axes
        from matplotlib.collections import LineCollection

        out = str(tmp_path / "chrono_circular_batched.png")
        svc = Chronogram(_make_args(plot_output=out, circular=True))
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        parent_map = build_parent_map(tree)
        root_to_tip = svc._compute_root_to_tip(tree)
        original_add_collection = matplotlib.axes.Axes.add_collection
        line_collections = []

        def count_collection(self, collection, *args, **kwargs):
            if isinstance(collection, LineCollection):
                line_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(
            matplotlib.axes.Axes, "add_collection", count_collection
        )

        svc._plot_circular(tree, parent_map, 1.0, root_to_tip)

        assert len(line_collections) >= 2
        assert any(
            len(collection.get_segments()) == 3
            and collection.get_linewidths()[0] == pytest.approx(0.8)
            for collection in line_collections
        )
        assert (tmp_path / "chrono_circular_batched.png").exists()

    def test_circular_timescale_reuses_boundary_circle_points(
        self, monkeypatch, tmp_path
    ):
        pytest.importorskip("matplotlib")
        import numpy as np

        out = str(tmp_path / "chrono_circular_timescale_points.png")
        svc = Chronogram(_make_args(plot_output=out, circular=True, root_age=20.0))
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        parent_map = build_parent_map(tree)
        root_to_tip = svc._compute_root_to_tip(tree)
        guide_calls = 0
        original_linspace = module.np.linspace

        def fake_timescale(_root_age, _timescale):
            return (
                [
                    ("ring1", 18.0, 15.0),
                    ("ring2", 15.0, 10.0),
                    ("ring3", 10.0, 5.0),
                    ("ring4", 5.0, 0.0),
                ],
                {},
            )

        def count_linspace(start, stop, num, *args, **kwargs):
            nonlocal guide_calls
            if start == 0 and stop == 2 * np.pi and num == 200:
                guide_calls += 1
            return original_linspace(start, stop, num, *args, **kwargs)

        monkeypatch.setattr(module, "get_timescale_for_range", fake_timescale)
        monkeypatch.setattr(module.np, "linspace", count_linspace)

        svc._plot_circular(tree, parent_map, 1.0, root_to_tip)

        assert guide_calls == 1
        assert (tmp_path / "chrono_circular_timescale_points.png").exists()

    def test_circular_plot_batches_hpd_intervals(self, monkeypatch, tmp_path):
        pytest.importorskip("matplotlib")
        import matplotlib.axes
        from matplotlib.collections import LineCollection

        out = str(tmp_path / "chrono_circular_hpd_batched.png")
        svc = Chronogram(_make_args(plot_output=out, circular=True))
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        parent_map = build_parent_map(tree)
        root_to_tip = svc._compute_root_to_tip(tree)
        hpd_intervals = {
            id(clade): (10.0, 20.0)
            for clade in svc._iter_preorder(tree.root)
            if clade.clades and clade != tree.root
        }
        original_plot = matplotlib.axes.Axes.plot
        original_add_collection = matplotlib.axes.Axes.add_collection
        line_collections = []

        def fail_hpd_plot(self, *args, **kwargs):
            if kwargs.get("color") == "#2b8cbe" and kwargs.get("lw") == 4:
                raise AssertionError("circular HPD intervals should be batched")
            return original_plot(self, *args, **kwargs)

        def count_collection(self, collection, *args, **kwargs):
            if isinstance(collection, LineCollection):
                line_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "plot", fail_hpd_plot)
        monkeypatch.setattr(
            matplotlib.axes.Axes, "add_collection", count_collection
        )

        svc._plot_circular(tree, parent_map, 1.0, root_to_tip, hpd_intervals)

        assert any(
            len(collection.get_segments()) == len(hpd_intervals)
            and collection.get_linewidths()[0] == 4
            for collection in line_collections
        )
        assert (tmp_path / "chrono_circular_hpd_batched.png").exists()

    def test_rectangular_creates_file(self, tmp_path):
        out = str(tmp_path / "chrono.png")
        svc = Chronogram(_make_args(plot_output=out))
        svc.run()
        assert (tmp_path / "chrono.png").exists()

    def test_circular_creates_file(self, tmp_path):
        out = str(tmp_path / "chrono_circ.png")
        svc = Chronogram(_make_args(plot_output=out, circular=True))
        svc.run()
        assert (tmp_path / "chrono_circ.png").exists()

    def test_node_ages_flag(self, tmp_path):
        out = str(tmp_path / "chrono_ages.png")
        svc = Chronogram(_make_args(plot_output=out, node_ages=True))
        svc.run()
        assert (tmp_path / "chrono_ages.png").exists()

    def test_timescale_epoch(self, tmp_path):
        out = str(tmp_path / "chrono_epoch.png")
        svc = Chronogram(_make_args(plot_output=out, timescale="epoch"))
        svc.run()
        assert (tmp_path / "chrono_epoch.png").exists()

    def test_timescale_period(self, tmp_path):
        out = str(tmp_path / "chrono_period.png")
        svc = Chronogram(_make_args(
            plot_output=out, timescale="period", root_age=250.0,
        ))
        svc.run()
        assert (tmp_path / "chrono_period.png").exists()

    def test_timescale_era(self, tmp_path):
        out = str(tmp_path / "chrono_era.png")
        svc = Chronogram(_make_args(
            plot_output=out, timescale="era", root_age=500.0,
        ))
        svc.run()
        assert (tmp_path / "chrono_era.png").exists()

    def test_json_output(self, tmp_path, capsys):
        out = str(tmp_path / "chrono.png")
        svc = Chronogram(_make_args(plot_output=out, json=True))
        svc.run()
        captured = capsys.readouterr()
        payload = json.loads(captured.out)
        assert payload["root_age"] == 70.0
        assert payload["n_tips"] == 8
        assert len(payload["node_ages"]) > 0
        # Check root node age is ~70
        root_entry = [n for n in payload["node_ages"] if n["n_descendants"] == 8]
        assert len(root_entry) == 1
        assert root_entry[0]["age_ma"] == pytest.approx(70.0, abs=0.1)

    def test_compute_root_to_tip_does_not_call_get_path(self, monkeypatch):
        svc = Chronogram(_make_args())
        tree = Phylo.read(StringIO("((A:1,B:2,D:3):3,C:4):0;"), "newick")

        def fail_get_path(self, *args, **kwargs):
            raise AssertionError("get_path should not be called")

        monkeypatch.setattr(type(tree), "get_path", fail_get_path)

        distances = svc._compute_root_to_tip(tree)

        by_name = {clade.name: distances[id(clade)] for clade in tree.get_terminals()}
        assert distances[id(tree.root)] == 0.0
        assert by_name == {"A": 4.0, "B": 5.0, "D": 6.0, "C": 4.0}

    def test_postorder_clades_direct_matches_biopython_order(self):
        tree = Phylo.read(
            StringIO("((A:1,B:2,C:3):4,(D:5,E:6):7,F:8):0;"),
            "newick",
        )

        direct = Chronogram._postorder_clades_direct(tree)
        reference = list(tree.find_clades(order="postorder"))

        assert [id(clade) for clade in direct] == [
            id(clade) for clade in reference
        ]

    def test_iter_preorder_preserves_order_without_reversed(self):
        class NoReversedList(list):
            def __reversed__(self):
                raise AssertionError("_iter_preorder should push children directly")

        root = Clade(name="root")
        left = Clade(name="left")
        middle = Clade(name="middle")
        right = Clade(name="right")
        left.clades = NoReversedList([Clade(name="left_a"), Clade(name="left_b")])
        middle.clades = NoReversedList([Clade(name="middle_a")])
        right.clades = NoReversedList()
        root.clades = NoReversedList([left, middle, right])

        order = [clade.name for clade in Chronogram._iter_preorder(root)]

        assert order == [
            "root",
            "left",
            "left_a",
            "left_b",
            "middle",
            "middle_a",
            "right",
        ]

    def test_preorder_clades_direct_preserves_order_without_reversed(self):
        class NoReversedList(list):
            def __reversed__(self):
                raise AssertionError(
                    "_preorder_clades_direct should push children directly"
                )

        root = Clade(name="root")
        left = Clade(name="left")
        middle = Clade(name="middle")
        right = Clade(name="right")
        left.clades = NoReversedList([Clade(name="left_a"), Clade(name="left_b")])
        middle.clades = NoReversedList([Clade(name="middle_a")])
        right.clades = NoReversedList()
        root.clades = NoReversedList([left, middle, right])

        tree = Tree(root=root)
        order = [clade.name for clade in Chronogram._preorder_clades_direct(tree)]

        assert order == [
            "root",
            "left",
            "left_a",
            "left_b",
            "middle",
            "middle_a",
            "right",
        ]

    def test_print_json_does_not_rewalk_descendant_terminals(
        self, monkeypatch
    ):
        svc = Chronogram(_make_args(json=True))
        tree = Phylo.read(StringIO("((A:1,B:2):3,C:4):0;"), "newick")
        root_to_tip = svc._compute_root_to_tip(tree)
        captured = {}

        def capture_payload(payload, **kwargs):
            captured["payload"] = payload

        def fail_get_terminals(self, *args, **kwargs):
            raise AssertionError("get_terminals should not be called")

        def fail_find_clades(self, *args, **kwargs):
            raise AssertionError("find_clades should not be called")

        monkeypatch.setattr(
            "phykit.services.tree.chronogram.print_json",
            capture_payload,
        )
        monkeypatch.setattr(type(tree.root), "get_terminals", fail_get_terminals)
        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        svc._print_json(tree, {}, 1.0, root_to_tip)

        payload = captured["payload"]
        assert payload["n_tips"] == 3
        root_entries = [
            node for node in payload["node_ages"]
            if node["n_descendants"] == 3
        ]
        assert len(root_entries) == 1
        assert root_entries[0]["descendants"] == ["A", "B", "C"]

    def test_ladderize(self, tmp_path):
        out = str(tmp_path / "chrono_ladder.png")
        svc = Chronogram(_make_args(plot_output=out, ladderize=True))
        svc.run()
        assert (tmp_path / "chrono_ladder.png").exists()

    def test_non_ultrametric_tree(self, tmp_path):
        """Non-ultrametric tree should still work (scaled by root age)."""
        out = str(tmp_path / "chrono_nonultra.png")
        svc = Chronogram(_make_args(
            tree="tests/sample_files/tree_simple.tre",
            plot_output=out, root_age=100.0,
        ))
        svc.run()
        assert (tmp_path / "chrono_nonultra.png").exists()

    def test_beast_hpd_auto_detected(self, tmp_path):
        """BEAST-style HPD annotations should be auto-detected and plotted."""
        out = str(tmp_path / "chrono_beast.png")
        svc = Chronogram(_make_args(
            tree="tests/sample_files/beast_annotated_tree.tre",
            plot_output=out, root_age=70.0,
        ))
        svc.run()
        assert (tmp_path / "chrono_beast.png").exists()

    def test_hpd_in_json(self, tmp_path, capsys):
        """JSON should include HPD intervals when present."""
        out = str(tmp_path / "chrono_beast.png")
        svc = Chronogram(_make_args(
            tree="tests/sample_files/beast_annotated_tree.tre",
            plot_output=out, root_age=70.0, json=True,
        ))
        svc.run()
        captured = capsys.readouterr()
        payload = json.loads(captured.out)
        # At least some nodes should have HPD
        hpd_nodes = [n for n in payload["node_ages"] if "hpd_lower" in n]
        assert len(hpd_nodes) > 0
        for n in hpd_nodes:
            assert n["hpd_lower"] < n["hpd_upper"]

    def test_beast_hpd_circular(self, tmp_path):
        """BEAST HPD should work in circular mode too."""
        out = str(tmp_path / "chrono_beast_circ.png")
        svc = Chronogram(_make_args(
            tree="tests/sample_files/beast_annotated_tree.tre",
            plot_output=out, root_age=70.0, circular=True,
        ))
        svc.run()
        assert (tmp_path / "chrono_beast_circ.png").exists()
