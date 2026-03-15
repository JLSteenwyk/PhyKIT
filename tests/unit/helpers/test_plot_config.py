import pytest

from phykit.helpers.plot_config import PlotConfig


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
