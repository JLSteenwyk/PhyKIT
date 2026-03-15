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
