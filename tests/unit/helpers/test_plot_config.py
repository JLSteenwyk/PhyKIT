import argparse

import pytest

from phykit.helpers.plot_config import PlotConfig, add_plot_arguments


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
