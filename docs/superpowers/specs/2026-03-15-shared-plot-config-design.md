# Shared Plot Configuration System

**Date:** 2026-03-15
**Status:** Approved
**Scope:** Phase 1 — `PlotConfig` + concatenation occupancy plot migration

## Problem

The concatenation occupancy plot (`--plot-occupancy` in `create_concatenation_matrix`) produces unreadable figures for large datasets. The figure height caps at 18 inches, font sizes are hardcoded at 7pt, and users have no way to customize visual parameters. This same problem exists across all ~26 plot-generating functions in PhyKIT.

A user requested: dynamic scaling based on dataset size, and CLI arguments to control title, legend position, dimensions, font sizes, and colors.

## Design

### New file: `phykit/helpers/plot_config.py`

#### `PlotConfig` dataclass

Fields (all optional except `dpi`, `show_title`):

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `fig_width` | `float \| None` | `None` (auto) | Figure width in inches |
| `fig_height` | `float \| None` | `None` (auto) | Figure height in inches |
| `dpi` | `int` | `300` | Output resolution |
| `show_title` | `bool` | `True` | Whether to display a title |
| `title` | `str \| None` | `None` | Custom title text (None = plot's default) |
| `legend_position` | `str \| None` | `None` | Matplotlib loc string, or `"none"` to hide |
| `label_fontsize` | `float \| None` | `None` (auto) | Font size for y-axis data labels (taxa names) |
| `xlabel_fontsize` | `float \| None` | `None` (auto) | Font size for x-axis data labels (gene names) |
| `title_fontsize` | `float \| None` | `None` (auto) | Font size for the title |
| `axis_fontsize` | `float \| None` | `None` (auto) | Font size for axis labels |
| `colors` | `list[str] \| None` | `None` | User-provided colors (hex or named) |

**Validation:** `fig_width`, `fig_height`, and `dpi` must be positive when provided. Fontsize fields must be non-negative: `0.0` is accepted as a special value meaning "hide these labels." `legend_position` is validated against matplotlib's accepted `loc` values plus the special `"none"` string. Invalid values produce a clear error message and exit.

**Flag precedence:** `--no-title` takes precedence over `--title`. If both are provided, the title is hidden.

#### `auto_scale(n_rows, n_cols)` classmethod

Returns a `PlotConfig` with computed defaults based on data dimensions.

When `n_rows` or `n_cols` is `None`, those dimensions are not used in scaling — the corresponding fields get static fallback defaults (e.g., `fig_height=8.0`, `fig_width=14.0`, `label_fontsize=7.0`).

- **fig_height:** `max(5.0, min(200.0, 3.0 + n_rows * 0.18))` — soft cap at 200 inches; a warning is printed if the uncapped value exceeds 200
- **fig_width:** `max(10.0, min(20.0, 8.0 + n_cols * 0.15))` when `n_cols is not None`; `14.0` otherwise
- **label_fontsize (y-axis):** `0.0` if `n_rows > 800` (meaning labels are hidden); otherwise `max(3.0, min(7.0, 7.0 - (n_rows - 50) * 0.008))` — scales from 7pt (<=50 rows) down to 3pt (~550+ rows). A `label_fontsize` of `0.0` signals the plot function to hide y-axis tick labels entirely.
- **xlabel_fontsize (x-axis):** `0.0` if `n_cols > 60` (labels hidden); otherwise `max(3.0, min(7.0, 7.0 - (n_cols - 20) * 0.1))` — scales from 7pt (<=20 cols) down to 3pt (~60 cols)
- **title_fontsize:** `12.0` (fixed default)
- **axis_fontsize:** `10.0` (fixed default)

These heuristics are starting points and will be tuned during implementation.

#### `from_args(args)` classmethod

Factory that reads CLI arguments from the argparse `Namespace` and returns a partially populated `PlotConfig` (user overrides only, auto-scaled fields left as `None`).

#### `resolve(n_rows=None, n_cols=None)` method

Called later when data dimensions are known. Fills any remaining `None` fields via `auto_scale(n_rows, n_cols)`. Returns `self` for chaining. Only fills fields that are still `None` — already-set values (from user CLI args or a prior `resolve` call) are not overwritten.

**Typical call flow:**
```python
# In __init__ / process_args (data dimensions not yet known):
self.plot_config = PlotConfig.from_args(args)

# Later, in the plot method (dimensions now known):
self.plot_config.resolve(n_rows=len(taxa), n_cols=len(genes))
```

#### `apply_to_figure(fig, ax, default_title, default_colors)` method

Applies config to a matplotlib figure. Responsibilities:
- Sets title (or hides it) using `title_fontsize`
- Positions or hides legend via `legend_position`
- Sets axis label font sizes using `axis_fontsize`
- Sets y-axis tick label font sizes using `label_fontsize` (hides labels if `0.0`)
- Sets x-axis tick label font sizes using `xlabel_fontsize` (hides labels if `0.0`)
- Returns resolved colors list (user-provided merged with `default_colors` fallback)

**Not** responsible for: `figsize` (set at figure creation time), `dpi` (set at `savefig` time), or plot-specific layout (gene boundaries, scatter points, etc.).

Each plot function calls this after drawing its content.

#### `add_plot_arguments(parser)` function

Registers the shared CLI argument group on any argparse parser:

```
--fig-width <float>       Figure width in inches (auto-scaled if omitted)
--fig-height <float>      Figure height in inches (auto-scaled if omitted)
--dpi <int>               Resolution in DPI (default: 300)
--no-title                Hide the plot title
--title <str>             Custom title text
--legend-position <str>   Legend location (e.g., "upper right", "lower left", "none")
--label-fontsize <float>  Font size for y-axis data labels (taxa names)
--xlabel-fontsize <float> Font size for x-axis data labels (gene names)
--title-fontsize <float>  Font size for the title
--axis-fontsize <float>   Font size for axis labels
--colors <str>            Comma-separated colors (hex or named)
```

**Color parsing:** The `--colors` string is split on commas in `from_args` and stored as `list[str]`. Whitespace around each color is stripped.

**Plot arg relevance:** These arguments are silently ignored when no plot is being generated (e.g., when `--plot-occupancy` is not set). This matches how `--plot-output` already behaves.

### Color handling

- Each plot function defines its own `default_colors` list.
- When the user provides `--colors`, those override from the left: color[0] replaces default[0], etc.
- If fewer user colors than defaults, remaining defaults are kept.
- If more user colors than defaults, extras are ignored.

Example for occupancy heatmap (3 colors: absent, gap, represented):
```
--colors "#000000,#cccccc,#e41a1c"
```

### Migration: concatenation occupancy plot

Changes to `create_concatenation_matrix.py`:
1. Import `PlotConfig` and `add_plot_arguments`
2. In `process_args`, call `PlotConfig.from_args(args)` and store as `self.plot_config`
3. In `_plot_concatenation_occupancy`, call `self.plot_config.resolve(n_rows=len(taxa), n_cols=len(alignment_paths))` then:
   - Replace `figsize=(14, fig_height)` with `(config.fig_width, config.fig_height)`
   - Replace hardcoded `ListedColormap` colors with resolved colors from `config.apply_to_figure()`
   - Replace hardcoded `fontsize=7` with `config.label_fontsize` / `config.xlabel_fontsize`
   - Use `config.apply_to_figure()` for title, legend, and tick label font sizes
   - Use `config.dpi` in `fig.savefig()`
   - Y-axis labels auto-hidden when `label_fontsize == 0.0` (handled by `apply_to_figure`)

Changes to `phykit.py`:
1. Call `add_plot_arguments(parser)` in the `create_concatenation_matrix` subcommand definition

### Rollout

**Phase 1 (this PR):**
- Create `phykit/helpers/plot_config.py`
- Migrate the concatenation occupancy plot
- Add shared CLI args to `create_concatenation_matrix` subcommand
- Unit tests for `PlotConfig`, `auto_scale`, `from_args`
- Update integration tests for the occupancy plot

**Phase 2 (follow-up PRs):**
- Migrate remaining ~25 plot functions in batches
- For `plot_alignment_qc`: keep existing `--width`/`--height`/`--dpi` as aliases for `--fig-width`/`--fig-height`/`--dpi` with a deprecation warning printed to stderr. Remove the aliases in a future major version.

## Files changed (Phase 1)

| File | Change |
|------|--------|
| `phykit/helpers/plot_config.py` | **New** — `PlotConfig`, `auto_scale`, `from_args`, `add_plot_arguments`, `apply_to_figure` |
| `phykit/services/alignment/create_concatenation_matrix.py` | Use `PlotConfig` in `_plot_concatenation_occupancy` and `process_args` |
| `phykit/phykit.py` | Add `add_plot_arguments()` call to `create_concatenation_matrix` parser |
| `tests/unit/helpers/test_plot_config.py` | **New** — unit tests for `PlotConfig` |
| `tests/unit/services/alignment/test_create_concatenation_matrix.py` | Update plot tests to cover new config |
