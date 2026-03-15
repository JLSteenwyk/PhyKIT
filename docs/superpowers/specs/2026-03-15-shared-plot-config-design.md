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
| `label_fontsize` | `float \| None` | `None` (auto) | Font size for data labels (taxa, genes) |
| `title_fontsize` | `float \| None` | `None` (auto) | Font size for the title |
| `axis_fontsize` | `float \| None` | `None` (auto) | Font size for axis labels |
| `colors` | `list[str] \| None` | `None` | User-provided colors (hex or named) |

#### `auto_scale(n_rows, n_cols)` classmethod

Returns a `PlotConfig` with computed defaults based on data dimensions:

- **fig_height:** `max(5.0, 3.0 + n_rows * 0.18)` — no upper cap, allowing large datasets to produce tall figures
- **fig_width:** `max(10.0, min(20.0, 8.0 + n_cols * 0.15))` when n_cols is relevant; `14.0` otherwise
- **label_fontsize:** `max(3.0, min(7.0, 7.0 - (n_rows - 50) * 0.008))` — scales from 7pt (<=50 rows) down to 3pt (~550+ rows); labels hidden beyond ~800 rows
- **title_fontsize:** `12.0` (fixed default)
- **axis_fontsize:** `10.0` (fixed default)

These heuristics are starting points and will be tuned during implementation.

#### `from_args(args, n_rows=None, n_cols=None)` classmethod

Factory that:
1. Reads CLI arguments from the argparse `Namespace`
2. For any field left as `None`, fills via `auto_scale(n_rows, n_cols)`
3. Returns a fully populated `PlotConfig`

#### `apply_to_figure(fig, ax, default_title, default_colors)` method

Applies config to a matplotlib figure:
- Sets title (or hides it) with `title_fontsize`
- Positions or hides legend via `legend_position`
- Returns resolved colors list (user-provided or `default_colors` fallback)

Each plot function calls this instead of hardcoding visual parameters.

#### `add_plot_arguments(parser)` function

Registers the shared CLI argument group on any argparse parser:

```
--fig-width <float>       Figure width in inches (auto-scaled if omitted)
--fig-height <float>      Figure height in inches (auto-scaled if omitted)
--dpi <int>               Resolution in DPI (default: 300)
--no-title                Hide the plot title
--title <str>             Custom title text
--legend-position <str>   Legend location (e.g., "upper right", "lower left", "none")
--label-fontsize <float>  Font size for data labels
--title-fontsize <float>  Font size for the title
--axis-fontsize <float>   Font size for axis labels
--colors <str>            Comma-separated colors (hex or named)
```

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
2. In `process_args`, build `PlotConfig` from args
3. In `_plot_concatenation_occupancy`:
   - Replace `figsize=(14, fig_height)` with `(config.fig_width, config.fig_height)`
   - Replace hardcoded `ListedColormap` colors with `config` colors
   - Replace hardcoded `fontsize=7` with `config.label_fontsize`
   - Use `config.apply_to_figure()` for title and legend
   - Use `config.dpi` in `fig.savefig()`
   - Auto-hide y-axis labels when taxa count exceeds threshold (~800)

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
- Deprecate `plot_alignment_qc`'s standalone `--width`/`--height`/`--dpi` in favor of shared args

## Files changed (Phase 1)

| File | Change |
|------|--------|
| `phykit/helpers/plot_config.py` | **New** — `PlotConfig`, `auto_scale`, `from_args`, `add_plot_arguments`, `apply_to_figure` |
| `phykit/services/alignment/create_concatenation_matrix.py` | Use `PlotConfig` in `_plot_concatenation_occupancy` and `process_args` |
| `phykit/phykit.py` | Add `add_plot_arguments()` call to `create_concatenation_matrix` parser |
| `tests/unit/helpers/test_plot_config.py` | **New** — unit tests for `PlotConfig` |
| `tests/unit/services/alignment/test_create_concatenation_matrix.py` | Update plot tests to cover new config |
