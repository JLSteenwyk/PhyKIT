# Circular Phylogram Layout Design Spec

## Overview

Add a `--circular` flag to all 11 phylogram-drawing commands in PhyKIT.
When set, the tree is drawn as a circular (radial/fan) phylogram with the
root at the center, branches radiating outward, and tips arranged around
the perimeter. Curved arcs connect sister clades at each internal node
(classic FigTree/iTOL style).

Combinable with `--cladogram` (equal branch lengths, all tips on outer
ring) and `--ladderize`.

`consensus_network` already uses a circular layout and ignores `--circular`.

## CLI Interface

The `--circular` flag is added to the shared `add_plot_arguments()` in
`PlotConfig`, making it available to all plot-generating commands:

```
phykit quartet_pie -t <tree> -g <genes> -o output.png --circular
phykit character_map -t <tree> -d <data> -o output.png --circular --cladogram
phykit phylo_heatmap -t <tree> -d <data> -o output.png --circular
```

Flag combinations:
- `--circular` alone = circular phylogram (radius = branch length)
- `--circular --cladogram` = circular cladogram (all tips at outer ring)
- `--circular --ladderize` = circular + ladderized ordering

## Architecture

### New File

`phykit/helpers/circular_layout.py` — shared circular drawing utilities.
All services call into this module; none implement circular drawing inline.

### Functions

```python
def compute_circular_coords(tree, node_x, parent_map):
    """Convert rectangular (node_x, tip_order) to circular coordinates.

    Args:
        tree: BioPython tree
        node_x: dict mapping node id -> x-coordinate (radius).
                From either phylogram (branch length) or cladogram (depth).
        parent_map: dict mapping node id -> parent clade

    Returns:
        Dict[node_id -> {"x": float, "y": float,
                          "angle": float, "radius": float}]

    Tip angles are evenly distributed around 2π (one per tip in tree
    traversal order). Internal node angles are computed as the midpoint
    of the angular range spanned by their descendant tips (not circular
    mean, to avoid the degenerate case when descendants span > 180°).
    The root's angle is irrelevant since its radius is 0.

    n_tips is derived from tree.count_terminals() — not passed as a
    parameter to avoid inconsistency.
    """

def draw_circular_branches(ax, tree, coords, parent_map,
                           color="#333", lw=1.5):
    """Draw radial + arc branches for the entire tree.

    For each internal node (which is the parent of its children):
    1. Draw an arc at the parent's radius connecting the angles of
       its children. This is the circular equivalent of the vertical
       connector in rectangular mode. The arc is always drawn via
       the shorter path; for subtrees spanning > 180°, the sweep
       goes the long way to stay on the correct side.
    2. For each child, draw a radial segment from the parent's radius
       to the child's radius at the child's angle (straight line
       along the radius).

    Arc direction: arcs sweep from the minimum child angle to the
    maximum child angle. When the angular span > π, the arc wraps
    around 0/2π by adding 2π to the end angle before sweeping.
    """

def draw_circular_tip_labels(ax, tree, coords, fontsize=9, offset=0.03):
    """Draw tip labels with angle-aware rotation and alignment.

    Labels on the right half (angle between -90° and 90°) are
    left-aligned and rotated to follow the angle. Labels on the
    left half are right-aligned and rotated 180° so text reads
    left-to-right.
    """

def draw_circular_colored_branch(ax, parent_coords, child_coords,
                                  color, lw=1.5):
    """Draw a single radial branch segment with a specific color.

    Used by services that color individual branches (discordance_asymmetry,
    rate_heterogeneity, etc.).
    """

def draw_circular_colored_arc(ax, center_x, center_y, radius,
                               start_angle, end_angle, color="#333", lw=1.5):
    """Draw a colored arc segment at a given radius between two angles.

    Used for coloring the arc portion of branches (the circular equivalent
    of the vertical connector). Services that color branches must call
    this for the arc AND draw_circular_colored_branch for the radial
    segment.

    Arc direction follows the same sweep rules as draw_circular_branches.
    """

def circular_branch_points(parent_coords, child_coords, n_points):
    """Return evenly spaced (x, y, angle) points along a radial segment.

    Used by character_map to position circles along branches, and by
    cont_map/density_map/stochastic_character_map to draw multi-segment
    colored branches.
    """

def radial_offset(angle, distance):
    """Return (dx, dy) offset for placing annotations along a radius.

    Used by character_map to position text labels radially outward
    (character index) and inward (state transition) from a point on
    a branch.

    Returns (cos(angle) * distance, sin(angle) * distance).
    """

def draw_circular_gradient_branch(ax, parent_coords, child_coords,
                                   colors_or_cmap, values, lw=1.5):
    """Draw a radial branch with a color gradient.

    Used by cont_map for continuous trait gradients along branches.
    Divides the radial segment into small sub-segments, each colored
    by interpolating the colormap.
    """

def draw_circular_multi_segment_branch(ax, parent_coords, child_coords,
                                        segments, state_colors, lw=1.5):
    """Draw a radial branch with multiple colored segments.

    Used by stochastic_character_map and density_map where branch
    history is a list of (fraction, state) segments.

    Args:
        segments: list of (start_frac, end_frac, state) along the branch
        state_colors: dict mapping state -> color
    """
```

### Coordinate System

The existing `node_x` values (computed by phylogram or cladogram mode)
are reinterpreted as radial distance:

```
radius = node_x[id]  (taken directly, scaled to fit figure)
angle = 2π * tip_index / n_tips  (for tips, evenly distributed)
angle = midpoint of descendant angular range  (for internal nodes)

x_plot = center_x + radius * cos(angle)
y_plot = center_y + radius * sin(angle)
```

Internal node angle computation uses the **midpoint of the angular range**
spanned by descendant tips, not the circular mean. This avoids the
degenerate case where `atan2(~0, ~0)` is undefined when descendants span
the full circle (e.g., the root node). For the root specifically, the
angle is irrelevant since radius = 0 maps to (center_x, center_y)
regardless.

The root's radius is taken from `node_x`, which is 0.0 in practice.
The code does not hard-code this assumption.

### Arc Drawing: Direction and Wrap-Around

Arcs are drawn at the parent's radius from child_min_angle to
child_max_angle. When the angular span between children exceeds π,
the arc must wrap around the 0/2π boundary. The rule:

```python
span = max_angle - min_angle
if span > π:
    # Children straddle the 0/2π boundary — sweep the long way
    # by going from max_angle to min_angle + 2π
    draw arc from max_angle to min_angle + 2π
else:
    draw arc from min_angle to max_angle
```

The arc is rendered as a dense polyline (60+ points) computed
parametrically: `(center + radius * cos(θ), center + radius * sin(θ))`
for θ swept between the start and end angles.

### Figure Dimensions for Circular Mode

When `--circular` is set and both `fig_width` and `fig_height` are
auto-scaled (not set by the user), `PlotConfig.resolve()` produces a
square figure:

```python
if self.circular and self.fig_width is None and self.fig_height is None:
    size = max(8.0, min(14.0, 3.0 + n_rows * 0.08))
    self.fig_width = size
    self.fig_height = size
```

This keeps the logic centralized in PlotConfig rather than requiring
each service to handle it.

### Modified Files

| File | Change |
|------|--------|
| `phykit/helpers/plot_config.py` | Add `circular: bool = False` to PlotConfig, `--circular` to `add_plot_arguments`, square auto-scaling in `resolve()` |
| `phykit/helpers/circular_layout.py` | New file with all shared circular drawing functions |
| `phykit/phykit.py` | Add `--circular` to all help text blocks |
| `docs/usage/index.rst` | Add `--circular` to option lists |

### Service Modifications (11 files)

Each service's plot method gets a circular branch:

```python
if self.plot_config.circular:
    from ...helpers.circular_layout import (
        compute_circular_coords,
        draw_circular_branches,
        draw_circular_tip_labels,
    )
    coords = compute_circular_coords(tree, node_x, parent_map)
    ax.set_aspect("equal")
    ax.axis("off")
    draw_circular_branches(ax, tree, coords, parent_map)
    draw_circular_tip_labels(ax, tree, coords)
    # ... service-specific annotations using coords ...
else:
    # ... existing rectangular drawing code ...
```

### Per-Service Annotation Details

**1. quartet_pie** — Pie chart insets at internal nodes.
Convert `coords[id]["x"], coords[id]["y"]` to figure-fraction coordinates
using `ax.transData.transform()` → `fig.transFigure.inverted().transform()`.
Place inset axes at those positions. No rotation needed (pies are circles).

**2. character_map** — Colored scatter points + text on branches.
Use `circular_branch_points()` to get evenly-spaced positions along each
radial segment. Place scatter markers at those points. Use `radial_offset()`
to compute `xytext` for `ax.annotate()`: character index offset radially
outward, state transition offset radially inward.

**3. cont_map** — Continuous trait gradient along branches.
Use `draw_circular_gradient_branch()` for each radial segment.
Use `draw_circular_colored_arc()` for the arc at the parent's radius,
colored by the parent node's trait value (same as the vertical connector
color in rectangular mode).

**4. density_map** — Posterior state density coloring on branches.
Use `draw_circular_multi_segment_branch()` with the density-averaged
state colors for radial segments.
Use `draw_circular_colored_arc()` for arcs, colored by parent state.

**5. stochastic_character_map** — Multi-segment colored branches.
Use `draw_circular_multi_segment_branch()` with the stochastic mapping
history segments for radial segments.
Use `draw_circular_colored_arc()` for arcs, colored by parent state.

**6. discordance_asymmetry** — Branches colored by asymmetry ratio.
Use `draw_circular_colored_branch()` per radial segment.
Use `draw_circular_colored_arc()` for arcs, both with the branch color.

**7. rate_heterogeneity** — Branches colored by regime.
Same as discordance_asymmetry — `draw_circular_colored_branch()` +
`draw_circular_colored_arc()`.

**8. ancestral_reconstruction** — Either contmap (gradient branches)
or discrete pie charts at nodes. Contmap uses `draw_circular_gradient_branch()`
+ `draw_circular_colored_arc()`.
Discrete pies use the same inset-axes approach as quartet_pie.

**9. concordance_asr** — Contmap-style gradient branches.
Uses `draw_circular_gradient_branch()` + `draw_circular_colored_arc()`.

**10. phylo_heatmap** — Tree + heatmap as concentric rings.
In circular mode, the tree is drawn in the inner portion and the heatmap
wraps around the outside as concentric rings (one ring per trait column).

Layout:
- Single axes with `ax.set_aspect("equal")` and `ax.axis("off")`
- Tree occupies the inner region (radius 0 to `inner_radius`)
- Each trait column becomes one concentric ring outside the tree
- `inner_radius` determined by `--split` (default 0.3 = 30% of total
  radius is tree, 70% is rings)
- Each taxon occupies an angular wedge of width `2π / n_tips`
- Wedges are drawn with `matplotlib.patches.Wedge(center, r_outer,
  theta1_deg, theta2_deg, width=ring_width)` colored by trait value
- Trait name labels positioned radially outside each ring
- Colorbar placed to the right of the figure via `fig.colorbar()`

**11. cophylo** — Two circular trees side by side.
In circular mode, draw two separate circular trees in a 1x2 subplot
grid (`fig, (ax1, ax2) = plt.subplots(1, 2)`). Each tree is drawn
as an independent circular phylogram. Connecting lines between matched
taxa are drawn using `matplotlib.patches.ConnectionPatch` which can
span across subplot boundaries. Lines connect from the tip position on
ax1 to the matched tip position on ax2.

### `consensus_network` Behavior

`consensus_network` already draws a circular splits graph using its own
polar coordinate system. The `--circular` flag is silently ignored by
this command (no error, no change in behavior).

### Matplotlib Setup for Circular Mode

```python
fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))
ax.set_aspect("equal")
ax.axis("off")  # no axes frame for circular layout
```

`PlotConfig.apply_to_figure()` detects `circular` mode and skips
axis label/tick adjustments (they are invisible with `axis("off")`):

```python
if self.circular:
    ax.set_aspect("equal")
    ax.axis("off")
    return self.merge_colors(default_colors)
```

## Text/JSON Output

No changes to text or JSON output — `--circular` is purely visual.
All computed values (CI, RI, quartet proportions, etc.) are identical
regardless of layout mode.

## Test Plan

### Unit Tests (`tests/unit/helpers/test_circular_layout.py`)

- `test_tip_angles_evenly_spaced` — n tips produce angles at 2π/n intervals
- `test_tip_angles_monotonic` — tip angles increase monotonically in traversal order
- `test_internal_node_angle_between_children` — internal angle is between min and max child angles
- `test_internal_angle_large_subtree` — correct angle for subtree spanning > 180°
- `test_root_at_center` — root coords map to (center_x, center_y) when radius=0
- `test_coords_cartesian_on_circle` — x² + y² ≈ radius² for each node
- `test_label_alignment_right_half` — text-anchor="start" for angles in [-90°, 90°]
- `test_label_alignment_left_half` — text-anchor="end" for angles outside [-90°, 90°]
- `test_branch_points_count` — circular_branch_points returns correct number
- `test_branch_points_on_radial_line` — all points lie on the radial path
- `test_radial_offset_direction` — offset points outward at correct angle
- `test_arc_wrap_around` — arc crossing 0/2π boundary drawn correctly
- `test_cladogram_circular_tips_same_radius` — with cladogram coords, all tips at max radius

### Integration Tests

Add one `test_*_circular` test to each existing integration test file
(11 total). Each test invokes the command with `--circular` and verifies:
- Output file exists and is non-empty
- No errors written to stderr

### Visual Validation

Generate sample circular figures for a few commands and visually inspect.
