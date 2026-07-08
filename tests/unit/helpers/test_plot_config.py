import argparse
from io import StringIO

import pytest
from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin

from phykit.helpers.plot_config import (
    PlotConfig,
    _assign_internal_y_from_children,
    _mean_child_y,
    _preorder_clades_direct,
    _terminal_clades_direct,
    add_plot_arguments,
    build_parent_map,
    compute_node_positions,
    draw_tree_branches,
    draw_tip_labels,
)


def test_module_import_does_not_import_typing():
    import subprocess
    import sys

    code = """
import sys
import phykit.helpers.plot_config as module
assert hasattr(module, "PlotConfig")
assert "typing" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_module_import_does_not_import_dataclasses():
    import subprocess
    import sys

    code = """
import sys
import phykit.helpers.plot_config as module
assert module.PlotConfig().dpi == 300
assert "dataclasses" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_assign_internal_y_binary_children_use_indexed_path():
    class IndexedOnlyList(list):
        def __iter__(self):
            raise AssertionError("binary child y setup should index children")

    left = type("Clade", (), {})()
    right = type("Clade", (), {})()
    parent = type("Clade", (), {})()
    parent.clades = IndexedOnlyList([left, right])
    node_y = {id(left): 2.0, id(right): 6.0}

    assert _assign_internal_y_from_children(parent, node_y) == pytest.approx(4.0)


class TestPlotConfigDefaults:
    def test_default_values(self):
        config = PlotConfig()
        assert config.fig_width is None
        assert config.fig_height is None
        assert config.dpi == 300
        assert config.show_title is True
        assert config.title is None
        assert config.legend_position is None
        assert config.ylabel_fontsize is None
        assert config.xlabel_fontsize is None
        assert config.title_fontsize is None
        assert config.axis_fontsize is None
        assert config.colors is None

    def test_custom_values(self):
        config = PlotConfig(fig_width=10.0, fig_height=8.0, dpi=150, colors=["red", "blue"])
        assert config.fig_width == 10.0
        assert config.fig_height == 8.0
        assert config.dpi == 150
        assert config.colors == ["red", "blue"]


class TestPlotConfigValidation:
    def test_validate_positive_fig_width(self):
        config = PlotConfig(fig_width=-5.0)
        with pytest.raises(SystemExit):
            config.validate()

    def test_validate_positive_fig_height(self):
        config = PlotConfig(fig_height=0.0)
        with pytest.raises(SystemExit):
            config.validate()

    def test_validate_positive_dpi(self):
        config = PlotConfig(dpi=0)
        with pytest.raises(SystemExit):
            config.validate()

    def test_validate_nonnegative_ylabel_fontsize(self):
        config = PlotConfig(ylabel_fontsize=-1.0)
        with pytest.raises(SystemExit):
            config.validate()

    def test_validate_zero_ylabel_fontsize_is_valid(self):
        config = PlotConfig(ylabel_fontsize=0.0)
        config.validate()  # should not raise

    def test_validate_nonnegative_xlabel_fontsize(self):
        config = PlotConfig(xlabel_fontsize=-1.0)
        with pytest.raises(SystemExit):
            config.validate()

    def test_validate_zero_xlabel_fontsize_is_valid(self):
        config = PlotConfig(xlabel_fontsize=0.0)
        config.validate()  # should not raise

    def test_validate_valid_legend_position(self):
        config = PlotConfig(legend_position="upper right")
        config.validate()  # should not raise

    def test_validate_none_legend_position(self):
        config = PlotConfig(legend_position="none")
        config.validate()  # should not raise

    def test_validate_invalid_legend_position(self):
        config = PlotConfig(legend_position="top middle")
        with pytest.raises(SystemExit):
            config.validate()

    def test_validate_none_fields_pass(self):
        config = PlotConfig()
        config.validate()  # should not raise


class TestAutoScale:
    def test_none_dimensions_returns_static_defaults(self):
        config = PlotConfig.auto_scale(n_rows=None, n_cols=None)
        assert config.fig_height == 8.0
        assert config.fig_width == 14.0
        assert config.ylabel_fontsize == 7.0
        assert config.xlabel_fontsize == 7.0
        assert config.title_fontsize == 12.0
        assert config.axis_fontsize == 10.0

    def test_small_dataset(self):
        config = PlotConfig.auto_scale(n_rows=20, n_cols=10)
        assert config.fig_height == pytest.approx(max(5.0, min(200.0, 3.0 + 20 * 0.18)))
        assert config.ylabel_fontsize == 7.0  # <=50 rows stays at 7
        assert config.xlabel_fontsize == 7.0  # <=20 cols stays at 7

    def test_medium_dataset_scales_fonts(self):
        config = PlotConfig.auto_scale(n_rows=300, n_cols=40)
        assert config.ylabel_fontsize < 7.0
        assert config.ylabel_fontsize >= 3.0
        assert config.xlabel_fontsize < 7.0
        assert config.xlabel_fontsize >= 3.0

    def test_large_dataset_hides_ylabel(self):
        config = PlotConfig.auto_scale(n_rows=900, n_cols=10)
        assert config.ylabel_fontsize == 0.0

    def test_many_cols_hides_xlabel(self):
        config = PlotConfig.auto_scale(n_rows=10, n_cols=70)
        assert config.xlabel_fontsize == 0.0

    def test_fig_height_soft_cap(self, capsys):
        config = PlotConfig.auto_scale(n_rows=2000, n_cols=None)
        assert config.fig_height == 200.0
        captured = capsys.readouterr()
        assert "exceeds" in captured.err.lower() or "warning" in captured.err.lower()

    def test_fig_width_scales_with_cols(self):
        config = PlotConfig.auto_scale(n_rows=10, n_cols=50)
        expected = max(10.0, min(20.0, 8.0 + 50 * 0.15))
        assert config.fig_width == pytest.approx(expected)

    def test_fig_width_none_cols_uses_default(self):
        config = PlotConfig.auto_scale(n_rows=10, n_cols=None)
        assert config.fig_width == 14.0


class TestResolve:
    def test_resolve_fills_none_fields(self):
        config = PlotConfig()
        config.resolve(n_rows=50, n_cols=20)
        assert config.fig_width is not None
        assert config.fig_height is not None
        assert config.ylabel_fontsize is not None
        assert config.xlabel_fontsize is not None
        assert config.title_fontsize is not None
        assert config.axis_fontsize is not None

    def test_resolve_preserves_user_overrides(self):
        config = PlotConfig(fig_width=20.0, ylabel_fontsize=5.0)
        config.resolve(n_rows=900, n_cols=100)
        assert config.fig_width == 20.0  # user override preserved
        assert config.ylabel_fontsize == 5.0  # user override preserved, not set to 0.0

    def test_resolve_returns_self(self):
        config = PlotConfig()
        result = config.resolve(n_rows=10, n_cols=10)
        assert result is config

    def test_resolve_does_not_overwrite_prior_resolve(self):
        config = PlotConfig()
        config.resolve(n_rows=10, n_cols=10)
        first_height = config.fig_height
        config.resolve(n_rows=500, n_cols=500)
        assert config.fig_height == first_height  # not overwritten


class TestAddPlotArguments:
    def _make_parser(self):
        parser = argparse.ArgumentParser()
        add_plot_arguments(parser)
        return parser

    def test_defaults_all_none(self):
        parser = self._make_parser()
        args = parser.parse_args([])
        config = PlotConfig.from_args(args)
        assert config.fig_width is None
        assert config.fig_height is None
        assert config.dpi == 300
        assert config.show_title is True
        assert config.title is None
        assert config.legend_position is None
        assert config.ylabel_fontsize is None
        assert config.xlabel_fontsize is None
        assert config.colors is None

    def test_all_args_provided(self):
        parser = self._make_parser()
        args = parser.parse_args([
            "--fig-width", "12",
            "--fig-height", "8",
            "--dpi", "150",
            "--no-title",
            "--title", "My Title",
            "--legend-position", "lower left",
            "--ylabel-fontsize", "5",
            "--xlabel-fontsize", "4",
            "--title-fontsize", "14",
            "--axis-fontsize", "11",
            "--colors", "#ff0000, #00ff00, #0000ff",
        ])
        config = PlotConfig.from_args(args)
        assert config.fig_width == 12.0
        assert config.fig_height == 8.0
        assert config.dpi == 150
        assert config.show_title is False  # --no-title wins
        assert config.title == "My Title"
        assert config.legend_position == "lower left"
        assert config.ylabel_fontsize == 5.0
        assert config.xlabel_fontsize == 4.0
        assert config.title_fontsize == 14.0
        assert config.axis_fontsize == 11.0
        assert config.colors == ["#ff0000", "#00ff00", "#0000ff"]

    def test_colors_empty_entries_preserved(self):
        parser = self._make_parser()
        args = parser.parse_args(["--colors", ",,#e41a1c"])
        config = PlotConfig.from_args(args)
        assert config.colors == ["", "", "#e41a1c"]

    def test_zero_fontsize_accepted(self):
        parser = self._make_parser()
        args = parser.parse_args(["--ylabel-fontsize", "0", "--xlabel-fontsize", "0"])
        config = PlotConfig.from_args(args)
        assert config.ylabel_fontsize == 0.0
        assert config.xlabel_fontsize == 0.0

    def test_from_args_missing_plot_attrs_uses_defaults(self):
        args = argparse.Namespace(alignment_list="x", prefix="y")
        config = PlotConfig.from_args(args)
        assert config.fig_width is None
        assert config.dpi == 300

    def test_from_args_rejects_negative_fig_width(self):
        parser = self._make_parser()
        args = parser.parse_args(["--fig-width", "-5"])
        with pytest.raises(SystemExit):
            PlotConfig.from_args(args)

    def test_from_args_rejects_negative_fontsize(self):
        parser = self._make_parser()
        args = parser.parse_args(["--ylabel-fontsize", "-1"])
        with pytest.raises(SystemExit):
            PlotConfig.from_args(args)

    def test_from_args_rejects_invalid_legend_position(self):
        parser = self._make_parser()
        args = parser.parse_args(["--legend-position", "top middle"])
        with pytest.raises(SystemExit):
            PlotConfig.from_args(args)


class TestMergeColors:
    def test_no_user_colors_returns_defaults(self):
        config = PlotConfig()
        result = config.merge_colors(["A", "B", "C"])
        assert result == ["A", "B", "C"]

    def test_full_override(self):
        config = PlotConfig(colors=["X", "Y", "Z"])
        result = config.merge_colors(["A", "B", "C"])
        assert result == ["X", "Y", "Z"]

    def test_partial_override(self):
        config = PlotConfig(colors=["X", "Y"])
        result = config.merge_colors(["A", "B", "C"])
        assert result == ["X", "Y", "C"]

    def test_extra_user_colors_ignored(self):
        config = PlotConfig(colors=["X", "Y", "Z", "W"])
        result = config.merge_colors(["A", "B", "C"])
        assert result == ["X", "Y", "Z"]

    def test_empty_entries_preserve_defaults(self):
        config = PlotConfig(colors=["", "", "#e41a1c"])
        result = config.merge_colors(["#525252", "#d9d9d9", "#2b8cbe"])
        assert result == ["#525252", "#d9d9d9", "#e41a1c"]

    def test_single_color_plot(self):
        config = PlotConfig(colors=["red"])
        result = config.merge_colors(["#2b8cbe"])
        assert result == ["red"]


class TestComputeNodePositions:
    def _make_tree(self):
        return Phylo.read(StringIO("((A:1,B:2):3,C:4);"), "newick")

    def test_mean_child_y_scans_children_once(self):
        class SinglePassChildren:
            def __init__(self, children):
                self.children = children
                self.iterations = 0

            def __iter__(self):
                self.iterations += 1
                yield from self.children

        children = [object(), object(), object()]
        iterable = SinglePassChildren(children)
        node_y = {
            id(children[0]): 1.0,
            id(children[2]): 5.0,
        }

        assert _mean_child_y(iterable, node_y) == pytest.approx(3.0)
        assert iterable.iterations == 1
        assert _mean_child_y([object()], {}) == 0.0

    def test_build_parent_map_uses_direct_tree_traversal(self, monkeypatch):
        tree = self._make_tree()

        def fail_find_clades(*args, **kwargs):
            raise AssertionError("generic tree traversal should not be used")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        parent_map = build_parent_map(tree)
        tips = [child for child in tree.root.clades[0].clades]

        assert parent_map[id(tree.root.clades[0])] is tree.root
        assert parent_map[id(tree.root.clades[1])] is tree.root
        assert parent_map[id(tips[0])] is tree.root.clades[0]
        assert parent_map[id(tips[1])] is tree.root.clades[0]

    def test_build_parent_map_handles_mixed_child_counts(self, monkeypatch):
        tree = Phylo.read(StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);"), "newick")

        def fail_find_clades(*args, **kwargs):
            raise AssertionError("generic tree traversal should not be used")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        parent_map = build_parent_map(tree)
        binary = tree.root.clades[1]
        trifurcation = tree.root.clades[2]

        assert parent_map[id(tree.root.clades[0])] is tree.root
        assert parent_map[id(binary)] is tree.root
        assert parent_map[id(trifurcation)] is tree.root
        assert parent_map[id(binary.clades[0])] is binary
        assert parent_map[id(binary.clades[1])] is binary
        assert parent_map[id(trifurcation.clades[0])] is trifurcation
        assert parent_map[id(trifurcation.clades[1])] is trifurcation
        assert parent_map[id(trifurcation.clades[2])] is trifurcation

    def test_compute_node_positions_phylogram_direct_traversal(self, monkeypatch):
        tree = self._make_tree()
        parent_map = build_parent_map(tree)
        tips = tree.get_terminals()
        internal = parent_map[id(tips[0])]

        def fail_traversal(*args, **kwargs):
            raise AssertionError("generic tree traversal should not be used")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)
        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)

        node_x, node_y = compute_node_positions(tree, parent_map)

        assert node_x[id(tree.root)] == pytest.approx(0.0)
        assert node_x[id(internal)] == pytest.approx(3.0)
        assert node_x[id(tips[0])] == pytest.approx(4.0)
        assert node_x[id(tips[1])] == pytest.approx(5.0)
        assert node_x[id(tips[2])] == pytest.approx(4.0)
        assert node_y[id(tips[0])] == 0
        assert node_y[id(tips[1])] == 1
        assert node_y[id(tips[2])] == 2
        assert node_y[id(internal)] == pytest.approx(0.5)
        assert node_y[id(tree.root)] == pytest.approx(1.25)

    def test_draw_tree_branches_uses_direct_traversal(self, monkeypatch):
        tree = self._make_tree()
        parent_map = build_parent_map(tree)
        node_x, node_y = compute_node_positions(tree, parent_map)

        def fail_find_clades(*args, **kwargs):
            raise AssertionError("generic tree traversal should not be used")

        class NoopAxes:
            def __init__(self):
                self.plot_calls = 0

            def plot(self, *args, **kwargs):
                self.plot_calls += 1

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        ax = NoopAxes()
        draw_tree_branches(ax, tree, node_x, node_y, parent_map)

        assert ax.plot_calls == 2 * len(parent_map)

    def test_draw_tree_branches_batches_real_matplotlib_axes(self, monkeypatch):
        pytest.importorskip("matplotlib")
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        tree = self._make_tree()
        parent_map = build_parent_map(tree)
        node_x, node_y = compute_node_positions(tree, parent_map)
        fig, ax = plt.subplots()

        def fail_plot(*args, **kwargs):
            raise AssertionError("real axes should batch branch drawing")

        monkeypatch.setattr(ax, "plot", fail_plot)

        try:
            draw_tree_branches(ax, tree, node_x, node_y, parent_map)
            assert len(ax.collections) == 2
        finally:
            plt.close(fig)

    def test_compute_node_positions_cladogram_direct_traversal(self, monkeypatch):
        tree = self._make_tree()
        parent_map = build_parent_map(tree)
        tips = tree.get_terminals()
        internal = parent_map[id(tips[0])]

        def fail_traversal(*args, **kwargs):
            raise AssertionError("generic tree traversal should not be used")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)
        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)

        node_x, node_y = compute_node_positions(tree, parent_map, cladogram=True)

        assert node_x[id(tree.root)] == pytest.approx(0.0)
        assert node_x[id(internal)] == pytest.approx(0.5)
        assert all(node_x[id(tip)] == pytest.approx(1.0) for tip in tips)
        assert node_y[id(tree.root)] == pytest.approx(1.25)

    @pytest.mark.parametrize("cladogram", [False, True])
    def test_compute_node_positions_reuses_precomputed_preorder(self, cladogram):
        tree = self._make_tree()
        parent_map = build_parent_map(tree)
        preorder_clades = list(tree.find_clades(order="preorder"))

        expected_x, expected_y = compute_node_positions(
            tree,
            parent_map,
            cladogram=cladogram,
        )
        observed_x, observed_y = compute_node_positions(
            tree,
            parent_map,
            cladogram=cladogram,
            preorder_clades=preorder_clades,
        )

        assert observed_x == pytest.approx(expected_x)
        assert observed_y == pytest.approx(expected_y)

    def test_direct_clade_helpers_preserve_mixed_child_order(self, monkeypatch):
        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1,G:1);"),
            "newick",
        )
        expected_preorder = list(tree.find_clades(order="preorder"))
        expected_terminals = tree.get_terminals()

        def fail_traversal(*args, **kwargs):
            raise AssertionError("generic tree traversal should not be used")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)

        assert _preorder_clades_direct(tree) == expected_preorder
        assert _terminal_clades_direct(tree) == expected_terminals

    def test_draw_tip_labels_uses_direct_terminal_traversal(self, monkeypatch):
        tree = self._make_tree()
        parent_map = build_parent_map(tree)
        node_x, node_y = compute_node_positions(tree, parent_map)

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("standard trees should not call get_terminals")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        class Ax:
            def text(self, *args, **kwargs):
                return None

        draw_tip_labels(Ax(), tree, node_x, node_y)


class TestApplyToFigure:
    @pytest.fixture(autouse=True)
    def _skip_no_matplotlib(self):
        pytest.importorskip("matplotlib")

    def _make_fig_ax(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.plot([0, 1], [0, 1])
        return fig, ax

    def test_sets_title(self):
        import matplotlib.pyplot as plt
        fig, ax = self._make_fig_ax()
        config = PlotConfig(show_title=True, title_fontsize=14.0, axis_fontsize=10.0,
                            ylabel_fontsize=7.0, xlabel_fontsize=7.0)
        config.apply_to_figure(fig, ax, default_title="Test Title", default_colors=["red"])
        assert ax.get_title() == "Test Title"
        plt.close(fig)

    def test_custom_title(self):
        import matplotlib.pyplot as plt
        fig, ax = self._make_fig_ax()
        config = PlotConfig(show_title=True, title="Custom", title_fontsize=14.0,
                            axis_fontsize=10.0, ylabel_fontsize=7.0, xlabel_fontsize=7.0)
        config.apply_to_figure(fig, ax, default_title="Default", default_colors=["red"])
        assert ax.get_title() == "Custom"
        plt.close(fig)

    def test_no_title(self):
        import matplotlib.pyplot as plt
        fig, ax = self._make_fig_ax()
        config = PlotConfig(show_title=False, title_fontsize=14.0, axis_fontsize=10.0,
                            ylabel_fontsize=7.0, xlabel_fontsize=7.0)
        config.apply_to_figure(fig, ax, default_title="Title", default_colors=["red"])
        assert ax.get_title() == ""
        plt.close(fig)

    def test_hides_legend_when_none(self):
        import matplotlib.pyplot as plt
        from matplotlib.patches import Patch
        fig, ax = self._make_fig_ax()
        ax.legend(handles=[Patch(facecolor="red", label="A")])
        config = PlotConfig(legend_position="none", title_fontsize=14.0, axis_fontsize=10.0,
                            ylabel_fontsize=7.0, xlabel_fontsize=7.0, show_title=True)
        config.apply_to_figure(fig, ax, default_title="T", default_colors=["red"])
        legend = ax.get_legend()
        assert legend is None or not legend.get_visible()
        plt.close(fig)

    def test_returns_merged_colors(self):
        import matplotlib.pyplot as plt
        fig, ax = self._make_fig_ax()
        config = PlotConfig(colors=["blue"], title_fontsize=14.0, axis_fontsize=10.0,
                            ylabel_fontsize=7.0, xlabel_fontsize=7.0, show_title=True)
        result = config.apply_to_figure(fig, ax, default_title="T", default_colors=["red", "green"])
        assert result == ["blue", "green"]
        plt.close(fig)

    def test_hides_ylabel_when_zero(self):
        import matplotlib.pyplot as plt
        fig, ax = self._make_fig_ax()
        ax.set_yticks([0, 1])
        ax.set_yticklabels(["A", "B"])
        config = PlotConfig(ylabel_fontsize=0.0, xlabel_fontsize=7.0, title_fontsize=14.0,
                            axis_fontsize=10.0, show_title=True)
        config.apply_to_figure(fig, ax, default_title="T", default_colors=["red"])
        assert all(t.get_text() == "" for t in ax.get_yticklabels())
        plt.close(fig)

    def test_hides_xlabel_when_zero(self):
        import matplotlib.pyplot as plt
        fig, ax = self._make_fig_ax()
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["G1", "G2"])
        config = PlotConfig(xlabel_fontsize=0.0, ylabel_fontsize=7.0, title_fontsize=14.0,
                            axis_fontsize=10.0, show_title=True)
        config.apply_to_figure(fig, ax, default_title="T", default_colors=["red"])
        assert all(t.get_text() == "" for t in ax.get_xticklabels())
        plt.close(fig)
