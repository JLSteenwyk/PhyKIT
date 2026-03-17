# Color File Annotation Design Spec

## Overview

Add a `--color-file <path>` shared plot argument to all phylogram-drawing
commands in PhyKIT. The color file is an iTOL-inspired TSV that lets users
color tip labels, highlight clades with transparent background bands, and
color clade branches — all without modifying the tree file.

Works with all existing commands and combines with `--circular`,
`--cladogram`, and `--ladderize`.

## CLI Interface

The `--color-file` argument is added to the shared `add_plot_arguments()`
in `PlotConfig`, making it available to all plot-generating commands:

```
phykit quartet_pie -t <tree> -g <genes> -o out.png --color-file clades.tsv
phykit character_map -t <tree> -d <data> -o out.png --color-file clades.tsv --circular
phykit phylo_heatmap -t <tree> -d <data> -o out.png --color-file clades.tsv
```

## Input Format

Tab-separated file with 3-4 columns. Lines starting with `#` are comments.

```
# node(s)            type     color      label (optional)
human                label    #e41a1c
chimp                label    #e41a1c
human,chimp,gorilla  range    #ffe0e0    Primates
cat,dog              range    #e0e0ff    Carnivora
bear,raccoon         clade    #0000ff
```

### Column Definitions

| Column | Description |
|--------|-------------|
| 1: node(s) | Single taxon name for `label` type. Comma-separated taxa for `range` and `clade` types (MRCA resolution). |
| 2: type | `label`, `range`, or `clade` |
| 3: color | Hex color (`#rrggbb`) or named color (`red`, `blue`, etc.) |
| 4: label | Optional legend label for `range` and `clade` entries. Shown in the figure legend. |

### Entry Types

**`label`** — Color a tip's label text.
- Field 1: single taxon name
- Applies to: all commands
- Effect: changes the tip label's font color

**`range`** — Colored background band behind a clade.
- Field 1: comma-separated taxa whose MRCA defines the clade
- Applies to: all commands (including trait-colored ones — drawn behind branches)
- Effect: transparent rectangle (rectangular mode) or wedge (circular mode) spanning the clade's extent, drawn at zorder below branches
- Alpha: 0.15 (subtle background, does not obscure branches or annotations)

**`clade`** — Color the branches of a clade.
- Field 1: comma-separated taxa whose MRCA defines the clade
- Applies to: commands that do NOT already color branches by trait values. Silently ignored by `cont_map`, `density_map`, `stochastic_character_map`, `discordance_asymmetry`, `rate_heterogeneity`, and `concordance_asr` (these commands color branches by computed values which take precedence).
- Effect: all branches within the clade (from MRCA to tips) drawn in the specified color

### MRCA Resolution

For `range` and `clade` types, the listed taxa define the clade via
their most recent common ancestor (MRCA). The user only needs to list
enough taxa to uniquely identify the MRCA — typically one taxon from
each "side" of the clade:

```
# These are equivalent if human and gorilla span the Primates clade:
human,gorilla           range    #ffe0e0    Primates
human,chimp,gorilla     range    #ffe0e0    Primates
```

If a listed taxon is not found in the tree, it is silently skipped.
If fewer than 2 valid taxa remain after filtering, the entry is skipped
with a warning to stderr.

## Architecture

### New File

`phykit/helpers/color_annotations.py` — parsing and MRCA resolution.

### Functions

```python
def parse_color_file(path: str) -> dict:
    """Parse a color annotation file.

    Returns:
        {
            "labels": {taxon_name: color, ...},
            "ranges": [(taxa_list, color, label_or_None), ...],
            "clades": [(taxa_list, color, label_or_None), ...],
        }
    """

def resolve_mrca(tree, taxa: list):
    """Find the MRCA clade of the given taxa in the tree.

    Returns the BioPython Clade object, or None if fewer than 2
    taxa are found in the tree.
    """

def get_clade_tip_ids(clade) -> set:
    """Return set of id(tip) for all terminal descendants of a clade."""

def apply_label_colors(ax, tree, coords_or_node_xy, label_colors, circular=False):
    """Recolor tip labels that have entries in label_colors.

    In practice, tip labels are drawn by each service or by
    draw_circular_tip_labels. This function finds the matplotlib
    Text artists on the axes and recolors them, OR the services
    pass label_colors to the tip-drawing code.
    """

def draw_range_rect(ax, tree, clade, color, node_x, node_y, alpha=0.15):
    """Draw a colored rectangle behind a clade in rectangular mode.

    The rectangle spans:
    - x: from the MRCA's node_x to slightly beyond max tip x
    - y: from the minimum to maximum node_y of the clade's tips
    Drawn at zorder=0 (behind branches).
    """

def draw_range_wedge(ax, tree, clade, color, coords, alpha=0.15):
    """Draw a colored wedge behind a clade in circular mode.

    The wedge spans the angular range of the clade's tips, from
    the MRCA's radius to slightly beyond the max tip radius.
    Uses matplotlib.patches.Wedge.
    Drawn at zorder=0 (behind branches).
    """

def get_clade_branch_ids(tree, clade, parent_map) -> set:
    """Return set of node ids for all branches within a clade.

    Includes the MRCA and all its descendants (used by 'clade'
    type to color branches).
    """
```

### Modified Files

| File | Change |
|------|--------|
| `phykit/helpers/plot_config.py` | Add `color_file: Optional[str] = None` to PlotConfig, `--color-file` to `add_plot_arguments`, read in `from_args` |
| `phykit/helpers/color_annotations.py` | New file with parsing and rendering functions |

### Service Modifications (11 files)

Each service's plot method adds color annotation support:

```python
# After drawing the tree but before service-specific annotations:
if self.plot_config.color_file:
    from ...helpers.color_annotations import (
        parse_color_file,
        resolve_mrca,
        draw_range_rect,  # or draw_range_wedge for circular
    )
    color_data = parse_color_file(self.plot_config.color_file)

    # Draw range backgrounds (behind everything)
    for taxa_list, color, label in color_data["ranges"]:
        mrca = resolve_mrca(tree, taxa_list)
        if mrca:
            if self.plot_config.circular:
                draw_range_wedge(ax, tree, mrca, color, coords)
            else:
                draw_range_rect(ax, tree, mrca, color, node_x, node_y)

    # Apply clade branch colors (only for non-trait-colored commands)
    # ... apply label colors to tip text ...
```

The `range` type is applied FIRST (zorder=0, behind branches). Then
branches are drawn normally. Then annotations. This ensures ranges
never obscure other visual elements.

### Commands That Ignore `clade` Type

These commands color branches by computed trait values. The `clade`
branch coloring is silently ignored (ranges and labels still apply):

- `cont_map`
- `density_map`
- `stochastic_character_map`
- `discordance_asymmetry`
- `rate_heterogeneity`
- `concordance_asr`

### Range Rendering Details

**Rectangular mode (`draw_range_rect`):**
```
x_min = node_x[id(mrca)]
x_max = max(node_x[id(tip)] for tip in mrca.get_terminals()) + padding
y_min = min(node_y[id(tip)] for tip in mrca.get_terminals()) - 0.4
y_max = max(node_y[id(tip)] for tip in mrca.get_terminals()) + 0.4

ax.add_patch(Rectangle(
    (x_min, y_min), x_max - x_min, y_max - y_min,
    facecolor=color, alpha=0.15, edgecolor="none", zorder=0,
))
```

**Circular mode (`draw_range_wedge`):**
```
tip_angles = [coords[id(tip)]["angle"] for tip in mrca.get_terminals()]
angle_min = min(tip_angles) - half_gap
angle_max = max(tip_angles) + half_gap
r_inner = coords[id(mrca)]["radius"]
r_outer = max(coords[id(tip)]["radius"] for tip in mrca.get_terminals()) + padding

ax.add_patch(Wedge(
    (0, 0), r_outer, degrees(angle_min), degrees(angle_max),
    width=r_outer - r_inner,
    facecolor=color, alpha=0.15, edgecolor="none", zorder=0,
))
```

### Legend for Ranges and Clades

If field 4 (label) is provided, add a legend entry:
```python
from matplotlib.patches import Patch
legend_handles.append(
    Patch(facecolor=color, alpha=0.3, label=label)
)
```

The range/clade legend is appended to any existing legend
(e.g., quartet_pie's concordance legend). Position follows
`--legend-position`.

### Label Coloring Implementation

Two approaches depending on when tip labels are drawn:

**Approach 1 (preferred): Pass colors to drawing code.**
Modify `draw_circular_tip_labels` to accept an optional
`label_colors: dict` parameter. Each service's rectangular
tip-label loop also checks for label colors.

**Approach 2: Post-hoc recolor.**
After all drawing, iterate over `ax.texts` and match text
content to taxon names, recoloring as needed. Fragile but
requires no changes to existing drawing code.

Use Approach 1: add `label_colors` parameter to the shared
circular drawing function and to each service's rectangular
tip label code.

## Text/JSON Output

No changes — `--color-file` is purely visual. If `--json` is used,
the color file is not reflected in the output.

## Test Plan

### Unit Tests (`tests/unit/helpers/test_color_annotations.py`)

- `test_parse_label_entry` — single taxon label parsed correctly
- `test_parse_range_entry` — comma-separated taxa + color + label
- `test_parse_clade_entry` — comma-separated taxa + color
- `test_parse_comments_skipped` — lines starting with # ignored
- `test_parse_empty_lines_skipped` — blank lines ignored
- `test_parse_missing_label_field` — 3-column entry (no label) parsed OK
- `test_resolve_mrca_two_taxa` — finds correct MRCA
- `test_resolve_mrca_missing_taxon` — skips unknown taxa gracefully
- `test_resolve_mrca_single_taxon` — returns None (need ≥2)
- `test_draw_range_rect_no_error` — smoke test on simple tree
- `test_draw_range_wedge_no_error` — smoke test on circular coords
- `test_get_clade_branch_ids` — returns correct node set

### Integration Tests

Add `test_*_color_file` to each service's integration test file
(11 total). Each test creates a small color file, invokes the command
with `--color-file`, and verifies the output file exists.

### Visual Validation

Generate sample figures with colored ranges for a few commands and
visually inspect.
