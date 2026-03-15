# Shared Plot Configuration Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create a shared `PlotConfig` system with auto-scaling and user-overridable CLI args, then migrate the concatenation occupancy plot to use it.

**Architecture:** A `PlotConfig` dataclass in `phykit/helpers/plot_config.py` holds all visual parameters. A `from_args()` classmethod reads CLI overrides, and `resolve()` fills remaining fields with auto-scaled values based on data dimensions. An `apply_to_figure()` method applies config to any matplotlib figure. An `add_plot_arguments()` function registers shared CLI args on any argparse parser.

**Tech Stack:** Python dataclasses, argparse, matplotlib

**Spec:** `docs/superpowers/specs/2026-03-15-shared-plot-config-design.md`

---

## Chunk 1: PlotConfig dataclass, auto_scale, validation

### Task 1: PlotConfig dataclass with defaults and validation

**Files:**
- Create: `tests/unit/helpers/test_plot_config.py`
- Create: `phykit/helpers/plot_config.py`

- [ ] **Step 1: Write tests for PlotConfig construction and validation**

```python
# tests/unit/helpers/test_plot_config.py
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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/unit/helpers/test_plot_config.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'phykit.helpers.plot_config'`

- [ ] **Step 3: Implement PlotConfig dataclass with validation**

```python
# phykit/helpers/plot_config.py
import sys
from dataclasses import dataclass, field
from typing import List, Optional


VALID_LEGEND_POSITIONS = frozenset([
    "best", "upper right", "upper left", "lower left", "lower right",
    "right", "center left", "center right", "lower center", "upper center",
    "center", "none",
])


@dataclass
class PlotConfig:
    fig_width: Optional[float] = None
    fig_height: Optional[float] = None
    dpi: int = 300
    show_title: bool = True
    title: Optional[str] = None
    legend_position: Optional[str] = None
    ylabel_fontsize: Optional[float] = None
    xlabel_fontsize: Optional[float] = None
    title_fontsize: Optional[float] = None
    axis_fontsize: Optional[float] = None
    colors: Optional[List[str]] = None

    def validate(self) -> None:
        if self.fig_width is not None and self.fig_width <= 0:
            print(f"Error: --fig-width must be positive, got {self.fig_width}")
            sys.exit(2)
        if self.fig_height is not None and self.fig_height <= 0:
            print(f"Error: --fig-height must be positive, got {self.fig_height}")
            sys.exit(2)
        if self.dpi <= 0:
            print(f"Error: --dpi must be positive, got {self.dpi}")
            sys.exit(2)
        for name, val in [
            ("--ylabel-fontsize", self.ylabel_fontsize),
            ("--xlabel-fontsize", self.xlabel_fontsize),
            ("--title-fontsize", self.title_fontsize),
            ("--axis-fontsize", self.axis_fontsize),
        ]:
            if val is not None and val < 0:
                print(f"Error: {name} must be non-negative, got {val}")
                sys.exit(2)
        if (
            self.legend_position is not None
            and self.legend_position not in VALID_LEGEND_POSITIONS
        ):
            print(
                f"Error: --legend-position '{self.legend_position}' is not valid. "
                f"Choose from: {', '.join(sorted(VALID_LEGEND_POSITIONS))}"
            )
            sys.exit(2)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/unit/helpers/test_plot_config.py -v`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add phykit/helpers/plot_config.py tests/unit/helpers/test_plot_config.py
git commit -m "feat: add PlotConfig dataclass with validation"
```

---

### Task 2: auto_scale classmethod

**Files:**
- Modify: `tests/unit/helpers/test_plot_config.py`
- Modify: `phykit/helpers/plot_config.py`

- [ ] **Step 1: Write tests for auto_scale**

Append to `tests/unit/helpers/test_plot_config.py`:

```python
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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/unit/helpers/test_plot_config.py::TestAutoScale -v`
Expected: FAIL — `AttributeError: type object 'PlotConfig' has no attribute 'auto_scale'`

- [ ] **Step 3: Implement auto_scale**

Add to `PlotConfig` class in `phykit/helpers/plot_config.py`:

```python
    @classmethod
    def auto_scale(cls, n_rows=None, n_cols=None) -> "PlotConfig":
        # fig_height
        if n_rows is not None:
            raw_height = 3.0 + n_rows * 0.18
            fig_height = max(5.0, min(200.0, raw_height))
            if raw_height > 200.0:
                print(
                    f"Warning: computed figure height ({raw_height:.0f} in) "
                    f"exceeds 200 in; capping at 200. Use --fig-height to override.",
                    file=sys.stderr,
                )
        else:
            fig_height = 8.0

        # fig_width
        if n_cols is not None:
            fig_width = max(10.0, min(20.0, 8.0 + n_cols * 0.15))
        else:
            fig_width = 14.0

        # ylabel_fontsize
        if n_rows is not None:
            if n_rows > 800:
                ylabel_fontsize = 0.0
            else:
                ylabel_fontsize = max(3.0, min(7.0, 7.0 - (n_rows - 50) * 0.008))
        else:
            ylabel_fontsize = 7.0

        # xlabel_fontsize
        if n_cols is not None:
            if n_cols > 60:
                xlabel_fontsize = 0.0
            else:
                xlabel_fontsize = max(3.0, min(7.0, 7.0 - (n_cols - 20) * 0.1))
        else:
            xlabel_fontsize = 7.0

        return cls(
            fig_width=fig_width,
            fig_height=fig_height,
            ylabel_fontsize=ylabel_fontsize,
            xlabel_fontsize=xlabel_fontsize,
            title_fontsize=12.0,
            axis_fontsize=10.0,
        )
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/unit/helpers/test_plot_config.py -v`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add phykit/helpers/plot_config.py tests/unit/helpers/test_plot_config.py
git commit -m "feat: add auto_scale classmethod to PlotConfig"
```

---

### Task 3: resolve method

**Files:**
- Modify: `tests/unit/helpers/test_plot_config.py`
- Modify: `phykit/helpers/plot_config.py`

- [ ] **Step 1: Write tests for resolve**

Append to `tests/unit/helpers/test_plot_config.py`:

```python
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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/unit/helpers/test_plot_config.py::TestResolve -v`
Expected: FAIL — `AttributeError: 'PlotConfig' object has no attribute 'resolve'`

- [ ] **Step 3: Implement resolve**

Add to `PlotConfig` class in `phykit/helpers/plot_config.py`:

```python
    def resolve(self, n_rows=None, n_cols=None) -> "PlotConfig":
        defaults = PlotConfig.auto_scale(n_rows, n_cols)
        for fld in [
            "fig_width", "fig_height", "ylabel_fontsize", "xlabel_fontsize",
            "title_fontsize", "axis_fontsize",
        ]:
            if getattr(self, fld) is None:
                setattr(self, fld, getattr(defaults, fld))
        return self
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/unit/helpers/test_plot_config.py -v`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add phykit/helpers/plot_config.py tests/unit/helpers/test_plot_config.py
git commit -m "feat: add resolve method to PlotConfig"
```

---

## Chunk 2: from_args, add_plot_arguments, color merging, apply_to_figure

### Task 4: add_plot_arguments and from_args

**Files:**
- Modify: `tests/unit/helpers/test_plot_config.py`
- Modify: `phykit/helpers/plot_config.py`

- [ ] **Step 1: Write tests for add_plot_arguments and from_args**

Add `import argparse` and `from phykit.helpers.plot_config import add_plot_arguments` to the **top of** `tests/unit/helpers/test_plot_config.py` (alongside the existing imports), then append the following test class:

```python
import argparse

from phykit.helpers.plot_config import add_plot_arguments


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
        """from_args works even if args namespace lacks plot fields (graceful fallback)."""
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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/unit/helpers/test_plot_config.py::TestAddPlotArguments -v`
Expected: FAIL — `ImportError: cannot import name 'add_plot_arguments'`

- [ ] **Step 3: Implement add_plot_arguments and from_args**

Add to `phykit/helpers/plot_config.py` (outside the class):

```python
def add_plot_arguments(parser) -> None:
    group = parser.add_argument_group("plot options")
    group.add_argument("--fig-width", type=float, default=None, help="Figure width in inches (auto-scaled if omitted)")
    group.add_argument("--fig-height", type=float, default=None, help="Figure height in inches (auto-scaled if omitted)")
    group.add_argument("--dpi", type=int, default=300, help="Resolution in DPI (default: 300)")
    group.add_argument("--no-title", action="store_true", default=False, help="Hide the plot title")
    group.add_argument("--title", type=str, default=None, help="Custom title text")
    group.add_argument("--legend-position", type=str, default=None, help="Legend location (e.g., 'upper right', 'none' to hide)")
    group.add_argument("--ylabel-fontsize", type=float, default=None, help="Font size for y-axis labels; 0 to hide")
    group.add_argument("--xlabel-fontsize", type=float, default=None, help="Font size for x-axis labels; 0 to hide")
    group.add_argument("--title-fontsize", type=float, default=None, help="Font size for the title")
    group.add_argument("--axis-fontsize", type=float, default=None, help="Font size for axis labels")
    group.add_argument("--colors", type=str, default=None, help="Comma-separated colors (hex or named)")
```

Add to `PlotConfig` class:

```python
    @classmethod
    def from_args(cls, args) -> "PlotConfig":
        colors_str = getattr(args, "colors", None)
        colors = None
        if colors_str is not None:
            colors = [c.strip() for c in colors_str.split(",")]

        no_title = getattr(args, "no_title", False)

        config = cls(
            fig_width=getattr(args, "fig_width", None),
            fig_height=getattr(args, "fig_height", None),
            dpi=getattr(args, "dpi", 300),
            show_title=not no_title,
            title=getattr(args, "title", None),
            legend_position=getattr(args, "legend_position", None),
            ylabel_fontsize=getattr(args, "ylabel_fontsize", None),
            xlabel_fontsize=getattr(args, "xlabel_fontsize", None),
            title_fontsize=getattr(args, "title_fontsize", None),
            axis_fontsize=getattr(args, "axis_fontsize", None),
            colors=colors,
        )
        config.validate()
        return config
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/unit/helpers/test_plot_config.py -v`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add phykit/helpers/plot_config.py tests/unit/helpers/test_plot_config.py
git commit -m "feat: add add_plot_arguments and from_args to PlotConfig"
```

---

### Task 5: Color merging logic

**Files:**
- Modify: `tests/unit/helpers/test_plot_config.py`
- Modify: `phykit/helpers/plot_config.py`

- [ ] **Step 1: Write tests for merge_colors**

Append to `tests/unit/helpers/test_plot_config.py`:

```python
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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/unit/helpers/test_plot_config.py::TestMergeColors -v`
Expected: FAIL — `AttributeError: 'PlotConfig' object has no attribute 'merge_colors'`

- [ ] **Step 3: Implement merge_colors**

Add to `PlotConfig` class in `phykit/helpers/plot_config.py`:

```python
    def merge_colors(self, defaults: List[str]) -> List[str]:
        if self.colors is None:
            return list(defaults)
        result = list(defaults)
        for i, user_color in enumerate(self.colors[:len(defaults)]):
            if user_color:  # non-empty string overrides
                result[i] = user_color
        return result
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/unit/helpers/test_plot_config.py -v`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add phykit/helpers/plot_config.py tests/unit/helpers/test_plot_config.py
git commit -m "feat: add merge_colors method to PlotConfig"
```

---

### Task 6: apply_to_figure method

**Files:**
- Modify: `tests/unit/helpers/test_plot_config.py`
- Modify: `phykit/helpers/plot_config.py`

- [ ] **Step 1: Write tests for apply_to_figure**

Append to `tests/unit/helpers/test_plot_config.py`:

```python
class TestApplyToFigure:
    @pytest.fixture(autouse=True)
    def _skip_no_matplotlib(self):
        pytest.importorskip("matplotlib")

    def _make_fig_ax(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        # Draw something so legend/title can be applied
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
        assert len(ax.get_yticklabels()) == 0 or all(
            t.get_text() == "" for t in ax.get_yticklabels()
        )
        plt.close(fig)

    def test_hides_xlabel_when_zero(self):
        import matplotlib.pyplot as plt
        fig, ax = self._make_fig_ax()
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["G1", "G2"])
        config = PlotConfig(xlabel_fontsize=0.0, ylabel_fontsize=7.0, title_fontsize=14.0,
                            axis_fontsize=10.0, show_title=True)
        config.apply_to_figure(fig, ax, default_title="T", default_colors=["red"])
        assert len(ax.get_xticklabels()) == 0 or all(
            t.get_text() == "" for t in ax.get_xticklabels()
        )
        plt.close(fig)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/unit/helpers/test_plot_config.py::TestApplyToFigure -v`
Expected: FAIL — `AttributeError: 'PlotConfig' object has no attribute 'apply_to_figure'`

- [ ] **Step 3: Implement apply_to_figure**

Add to `PlotConfig` class in `phykit/helpers/plot_config.py`:

```python
    def apply_to_figure(self, fig, ax, default_title, default_colors):
        # Title
        if self.show_title:
            title_text = self.title if self.title is not None else default_title
            ax.set_title(title_text, fontsize=self.title_fontsize)
        else:
            ax.set_title("")

        # Legend
        if self.legend_position is not None:
            if self.legend_position == "none":
                legend = ax.get_legend()
                if legend is not None:
                    legend.set_visible(False)
            else:
                legend = ax.get_legend()
                if legend is not None:
                    handles = legend.legend_handles
                    labels = [t.get_text() for t in legend.get_texts()]
                    legend.remove()
                    ax.legend(handles=handles, labels=labels, loc=self.legend_position)

        # Axis label font sizes
        if self.axis_fontsize is not None:
            ax.xaxis.label.set_fontsize(self.axis_fontsize)
            ax.yaxis.label.set_fontsize(self.axis_fontsize)

        # Y-axis tick labels
        if self.ylabel_fontsize is not None:
            if self.ylabel_fontsize == 0.0:
                ax.set_yticklabels([])
            else:
                for label in ax.get_yticklabels():
                    label.set_fontsize(self.ylabel_fontsize)

        # X-axis tick labels
        if self.xlabel_fontsize is not None:
            if self.xlabel_fontsize == 0.0:
                ax.set_xticklabels([])
            else:
                for label in ax.get_xticklabels():
                    label.set_fontsize(self.xlabel_fontsize)

        # Colors
        return self.merge_colors(default_colors)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/unit/helpers/test_plot_config.py -v`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add phykit/helpers/plot_config.py tests/unit/helpers/test_plot_config.py
git commit -m "feat: add apply_to_figure method to PlotConfig"
```

---

## Chunk 3: Migrate concatenation occupancy plot and CLI

### Task 7: Register shared plot arguments on create_concatenation_matrix CLI

**Files:**
- Modify: `phykit/phykit.py:5166-5171`

- [ ] **Step 1: Add add_plot_arguments call**

In `phykit/phykit.py`, add an import near the top of the file (with other helper imports, matching the existing relative import style):

```python
from .helpers.plot_config import add_plot_arguments
```

Then after line 5170 (`parser.add_argument("--plot-output", ...)`), add:

```python
        add_plot_arguments(parser)
```

- [ ] **Step 2: Update the usage text in the docstring**

In `phykit/phykit.py`, update the usage section of `create_concatenation_matrix` (around lines 5129-5162) to include the new plot options in the help text:

Replace:
```
                Usage:
                phykit create_concatenation_matrix -a <file> -p <string>
                  [--threshold <float>] [--plot-occupancy]
                  [--plot-output <path>] [--json]
```

With:
```
                Usage:
                phykit create_concatenation_matrix -a <file> -p <string>
                  [--threshold <float>] [--plot-occupancy]
                  [--plot-output <path>] [--json]
                  [--fig-width <float>] [--fig-height <float>]
                  [--dpi <int>] [--no-title] [--title <str>]
                  [--legend-position <str>]
                  [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                  [--title-fontsize <float>] [--axis-fontsize <float>]
                  [--colors <str>]
```

And add after the `--plot-output` description block (after line 5159):
```
                --fig-width                 figure width in inches
                                            (auto-scaled if omitted)

                --fig-height                figure height in inches
                                            (auto-scaled if omitted)

                --dpi                       resolution in DPI
                                            (default: 300)

                --no-title                  hide the plot title

                --title                     custom title text

                --legend-position           legend location (e.g.,
                                            "upper right", "none")

                --ylabel-fontsize           font size for y-axis labels;
                                            0 to hide

                --xlabel-fontsize           font size for x-axis labels;
                                            0 to hide

                --title-fontsize            font size for the title

                --axis-fontsize             font size for axis labels

                --colors                    comma-separated colors
                                            (hex or named, e.g.,
                                            "#ff0000,blue,#00ff00")
```

- [ ] **Step 3: Verify the CLI parses new args**

Run: `python -m phykit create_concatenation_matrix --help`
Expected: Output includes `--fig-width`, `--ylabel-fontsize`, etc.

- [ ] **Step 4: Commit**

```bash
git add phykit/phykit.py
git commit -m "feat: register shared plot arguments on create_concatenation_matrix CLI"
```

---

### Task 8: Wire PlotConfig into CreateConcatenationMatrix

**Files:**
- Modify: `phykit/services/alignment/create_concatenation_matrix.py:1-28, 36-44`
- Modify: `tests/unit/services/alignment/test_create_concatenation_matrix.py`

- [ ] **Step 1: Update existing tests to include plot_config**

In `tests/unit/services/alignment/test_create_concatenation_matrix.py`, update the `args` fixture at line 23:

```python
@pytest.fixture
def args():
    return Namespace(alignment_list="/some/path/to/file", prefix="some_prefix")
```

This stays the same — `from_args` uses `getattr` with defaults, so missing attrs are fine.

Update `test_init_sets_alignment_list_path` to also check `plot_config`:

```python
    def test_init_sets_alignment_list_path(self, args):
        ccm = CreateConcatenationMatrix(args)
        assert ccm.alignment_list_path == args.alignment_list
        assert ccm.prefix == args.prefix
        assert ccm.output_file_path is None
        assert ccm.json_output is False
        assert ccm.plot_occupancy is False
        assert ccm.plot_output is None
        assert ccm.threshold == 0
        assert ccm.plot_config is not None
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/unit/services/alignment/test_create_concatenation_matrix.py::TestCreateConcatenationMatrix::test_init_sets_alignment_list_path -v`
Expected: FAIL — `AssertionError: ... has no attribute 'plot_config'`

- [ ] **Step 3: Wire PlotConfig into the service class**

In `phykit/services/alignment/create_concatenation_matrix.py`, add import at top:

```python
from ...helpers.plot_config import PlotConfig
```

Update `__init__` (lines 19-28) to store plot_config:

```python
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            alignment_list_path=parsed["alignment_list_path"],
            prefix=parsed["prefix"],
        )
        self.json_output = parsed["json_output"]
        self.plot_occupancy = parsed["plot_occupancy"]
        self.plot_output = parsed["plot_output"]
        self.threshold = parsed["threshold"]
        self.plot_config = parsed["plot_config"]
```

Update `process_args` (lines 36-44) to include plot_config:

```python
    def process_args(self, args) -> Dict[str, str]:
        return dict(
            alignment_list_path=args.alignment_list,
            prefix=args.prefix,
            json_output=getattr(args, "json", False),
            plot_occupancy=getattr(args, "plot_occupancy", False),
            plot_output=getattr(args, "plot_output", None),
            threshold=getattr(args, "threshold", 0),
            plot_config=PlotConfig.from_args(args),
        )
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/unit/services/alignment/test_create_concatenation_matrix.py -v`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add phykit/services/alignment/create_concatenation_matrix.py tests/unit/services/alignment/test_create_concatenation_matrix.py
git commit -m "feat: wire PlotConfig into CreateConcatenationMatrix"
```

---

### Task 9: Migrate _plot_concatenation_occupancy to use PlotConfig

> **Note:** This migration intentionally changes x-label visibility behavior. The original code hid gene labels when `len(alignment_paths) > 40`. The new code uses `auto_scale` which hides labels when `n_cols > 60` (via `xlabel_fontsize=0.0`). This means datasets with 41-60 genes will now show labels at a scaled-down font size, which is an improvement.

**Files:**
- Modify: `phykit/services/alignment/create_concatenation_matrix.py:46-129`
- Modify: `tests/unit/services/alignment/test_create_concatenation_matrix.py`

- [ ] **Step 1: Update the existing plot test**

In `tests/unit/services/alignment/test_create_concatenation_matrix.py`, update `test_plot_concatenation_occupancy` (line 234):

```python
    def test_plot_concatenation_occupancy(self, tmp_path, args):
        pytest.importorskip("matplotlib")
        ccm = CreateConcatenationMatrix(args)
        output_file = tmp_path / "occ.png"
        ccm._plot_concatenation_occupancy(
            taxa=["A", "B"],
            alignment_paths=["g1.fa", "g2.fa"],
            concatenated_seqs={"A": ["AC", "GT"], "B": ["A-", "G?"]},
            present_taxa_by_gene=[{"A", "B"}, {"A", "B"}],
            gene_lengths=[2, 2],
            output_file=str(output_file),
        )
        assert output_file.exists()

    def test_plot_concatenation_occupancy_pdf_output(self, tmp_path, args):
        pytest.importorskip("matplotlib")
        ccm = CreateConcatenationMatrix(args)
        output_file = tmp_path / "occ.pdf"
        ccm._plot_concatenation_occupancy(
            taxa=["A", "B"],
            alignment_paths=["g1.fa", "g2.fa"],
            concatenated_seqs={"A": ["AC", "GT"], "B": ["A-", "G?"]},
            present_taxa_by_gene=[{"A", "B"}, {"A", "B"}],
            gene_lengths=[2, 2],
            output_file=str(output_file),
        )
        assert output_file.exists()

    def test_plot_concatenation_occupancy_custom_config(self, tmp_path):
        pytest.importorskip("matplotlib")
        from phykit.helpers.plot_config import PlotConfig
        custom_args = Namespace(
            alignment_list="/x", prefix="x",
            fig_width=10.0, fig_height=6.0, dpi=72, no_title=True,
            title=None, legend_position="lower left",
            ylabel_fontsize=5.0, xlabel_fontsize=4.0,
            title_fontsize=10.0, axis_fontsize=8.0,
            colors="#000000,#ffffff,#ff0000",
        )
        ccm = CreateConcatenationMatrix(custom_args)
        output_file = tmp_path / "custom.png"
        ccm._plot_concatenation_occupancy(
            taxa=["A", "B", "C"],
            alignment_paths=["g1.fa", "g2.fa"],
            concatenated_seqs={"A": ["AC", "GT"], "B": ["A-", "G?"], "C": ["TT", "AA"]},
            present_taxa_by_gene=[{"A", "B"}, {"A", "C"}],
            gene_lengths=[2, 2],
            output_file=str(output_file),
        )
        assert output_file.exists()

    def test_plot_concatenation_occupancy_svg_output(self, tmp_path, args):
        pytest.importorskip("matplotlib")
        ccm = CreateConcatenationMatrix(args)
        output_file = tmp_path / "occ.svg"
        ccm._plot_concatenation_occupancy(
            taxa=["A", "B"],
            alignment_paths=["g1.fa", "g2.fa"],
            concatenated_seqs={"A": ["AC", "GT"], "B": ["A-", "G?"]},
            present_taxa_by_gene=[{"A", "B"}, {"A", "B"}],
            gene_lengths=[2, 2],
            output_file=str(output_file),
        )
        assert output_file.exists()
```

- [ ] **Step 2: Run the new tests to see the baseline**

Run: `pytest tests/unit/services/alignment/test_create_concatenation_matrix.py::TestCreateConcatenationMatrix::test_plot_concatenation_occupancy_custom_config -v`
Expected: May PASS (the function works but ignores config) or FAIL (if config wiring causes issues). Either way, this establishes baseline.

- [ ] **Step 3: Rewrite _plot_concatenation_occupancy to use PlotConfig**

Replace the method at lines 46-129 of `phykit/services/alignment/create_concatenation_matrix.py`:

```python
    def _plot_concatenation_occupancy(
        self,
        taxa: List[str],
        alignment_paths: List[str],
        concatenated_seqs: Dict[str, List[str]],
        present_taxa_by_gene: List[set],
        gene_lengths: List[int],
        output_file: str,
    ) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.colors import ListedColormap
            from matplotlib.patches import Patch
        except ImportError:
            print("matplotlib is required for --plot-occupancy. Install matplotlib and retry.")
            sys.exit(2)

        config = self.plot_config
        config.resolve(n_rows=len(taxa), n_cols=len(alignment_paths))

        # 0: absent gene block, 1: present but gap/ambiguous character, 2: represented character
        total_len = int(sum(gene_lengths))
        state_matrix = np.zeros((len(taxa), total_len), dtype=np.uint8)
        invalid_chars = set(["-", "?", "*", "X", "x", "N", "n"])

        gene_boundaries = []
        cursor = 0
        for gene_idx, gene_len in enumerate(gene_lengths):
            start = cursor
            end = cursor + gene_len
            gene_boundaries.append(end)
            present_taxa = present_taxa_by_gene[gene_idx]

            for taxon_idx, taxon in enumerate(taxa):
                if taxon not in present_taxa:
                    state_matrix[taxon_idx, start:end] = 0
                    continue
                seq = concatenated_seqs[taxon][gene_idx]
                for pos_idx, char in enumerate(seq):
                    if char in invalid_chars:
                        state_matrix[taxon_idx, start + pos_idx] = 1
                    else:
                        state_matrix[taxon_idx, start + pos_idx] = 2
            cursor = end

        # Sort taxa by total represented occupancy (state == 2), descending
        represented_counts = np.sum(state_matrix == 2, axis=1)
        order = np.argsort(-represented_counts)
        state_matrix = state_matrix[order, :]
        taxa_sorted = [taxa[idx] for idx in order]

        default_colors = ["#525252", "#d9d9d9", "#2b8cbe"]
        colors = config.merge_colors(default_colors)

        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))
        cmap = ListedColormap(colors)
        ax.imshow(state_matrix, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=2)

        for boundary in gene_boundaries[:-1]:
            ax.axvline(boundary - 0.5, color="black", linewidth=0.6, alpha=0.8)

        # Label genes at centers when feasible (controlled by xlabel_fontsize)
        if config.xlabel_fontsize and config.xlabel_fontsize > 0:
            starts = [0] + gene_boundaries[:-1]
            centers = [((s + e) / 2) - 0.5 for s, e in zip(starts, gene_boundaries)]
            labels = [os.path.basename(path) for path in alignment_paths]
            ax.set_xticks(centers)
            ax.set_xticklabels(labels, rotation=90, fontsize=config.xlabel_fontsize)
        else:
            ax.set_xticks([])
            ax.set_xlabel("Concatenated alignment positions (gene boundaries shown)")

        # Y-axis tick labels (controlled by ylabel_fontsize)
        if config.ylabel_fontsize and config.ylabel_fontsize > 0:
            ax.set_yticks(np.arange(len(taxa_sorted)))
            ax.set_yticklabels(taxa_sorted, fontsize=config.ylabel_fontsize)
        else:
            ax.set_yticks([])

        ax.set_ylabel("Taxa (sorted by represented occupancy)")

        legend_handles = [
            Patch(facecolor=colors[2], label="Represented character"),
            Patch(facecolor=colors[1], label="Gap/Ambiguous in present gene"),
            Patch(facecolor=colors[0], label="Gene absent (placeholder block)"),
        ]

        legend_loc = config.legend_position or "upper right"
        if legend_loc != "none":
            ax.legend(handles=legend_handles, loc=legend_loc, fontsize=8, frameon=True)

        # Apply title via config
        if config.show_title:
            title_text = config.title if config.title is not None else "Concatenation Occupancy Map"
            ax.set_title(title_text, fontsize=config.title_fontsize)

        if config.axis_fontsize:
            ax.xaxis.label.set_fontsize(config.axis_fontsize)
            ax.yaxis.label.set_fontsize(config.axis_fontsize)

        fig.tight_layout()
        fig.savefig(output_file, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)
```

- [ ] **Step 4: Run all tests to verify they pass**

Run: `pytest tests/unit/services/alignment/test_create_concatenation_matrix.py -v`
Expected: All PASS

- [ ] **Step 5: Run the integration tests too**

Run: `pytest tests/integration/alignment/test_create_concatenation_matrix_integration.py -v`
Expected: All PASS

- [ ] **Step 6: Commit**

```bash
git add phykit/services/alignment/create_concatenation_matrix.py tests/unit/services/alignment/test_create_concatenation_matrix.py
git commit -m "feat: migrate concatenation occupancy plot to use PlotConfig"
```

---

### Task 10: Final verification

- [ ] **Step 1: Run the full test suite**

Run: `pytest tests/ -v --timeout=120`
Expected: All existing tests pass, no regressions.

- [ ] **Step 2: Smoke test the CLI end-to-end**

Create two small test FASTA files and run the full command:

```bash
mkdir -p /tmp/phykit_test
echo -e ">A\nACGT\n>B\nA-GT" > /tmp/phykit_test/g1.fa
echo -e ">A\nTT\n>C\nTA" > /tmp/phykit_test/g2.fa
echo -e "/tmp/phykit_test/g1.fa\n/tmp/phykit_test/g2.fa" > /tmp/phykit_test/alignments.txt

# Test with default auto-scaling
python -m phykit create_concatenation_matrix -a /tmp/phykit_test/alignments.txt -p /tmp/phykit_test/out --plot-occupancy

# Test with custom config and PDF output
python -m phykit create_concatenation_matrix -a /tmp/phykit_test/alignments.txt -p /tmp/phykit_test/out2 --plot-occupancy --plot-output /tmp/phykit_test/out2.pdf --fig-width 10 --fig-height 6 --no-title --ylabel-fontsize 5 --colors "#000000,#cccccc,#e41a1c"

# Test with SVG output
python -m phykit create_concatenation_matrix -a /tmp/phykit_test/alignments.txt -p /tmp/phykit_test/out3 --plot-occupancy --plot-output /tmp/phykit_test/out3.svg

# Verify files exist
ls -la /tmp/phykit_test/out.occupancy.png /tmp/phykit_test/out2.pdf /tmp/phykit_test/out3.svg
```

Expected: All three files created successfully.

- [ ] **Step 3: Clean up temp files**

```bash
rm -rf /tmp/phykit_test
```
