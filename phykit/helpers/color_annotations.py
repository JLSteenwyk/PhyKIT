"""Shared utilities for parsing color annotation files and drawing
colored ranges/clades on phylogenetic tree plots."""

import sys
from math import atan2, degrees, pi

from phykit.errors import PhykitUserError


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------


def parse_color_file(path: str) -> dict:
    """Read a TSV color-annotation file.

    Expected columns (tab-separated):
        field1  field2  field3  [field4]
        node    type    color   [label]

    *type* is one of ``label``, ``range``, or ``clade``.

    Returns a dict with keys ``labels``, ``ranges``, ``clades``.
    """
    try:
        with open(path) as fh:
            lines = fh.readlines()
    except FileNotFoundError:
        raise PhykitUserError([f"Color file not found: {path}"])

    labels: dict = {}
    ranges: list = []
    clades: list = []

    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        parts = stripped.split("\t")
        if len(parts) < 3:
            continue

        field1 = parts[0].strip()
        field2 = parts[1].strip().lower()
        field3 = parts[2].strip()
        field4 = parts[3].strip() if len(parts) > 3 else None

        if field2 == "label":
            labels[field1] = field3
        elif field2 == "range":
            taxa = [t.strip() for t in field1.split(",")]
            ranges.append((taxa, field3, field4))
        elif field2 == "clade":
            taxa = [t.strip() for t in field1.split(",")]
            clades.append((taxa, field3, field4))

    return {
        "labels": labels,
        "ranges": ranges,
        "clades": clades,
    }


# ---------------------------------------------------------------------------
# Tree helpers
# ---------------------------------------------------------------------------


def resolve_mrca(tree, taxa: list):
    """Return the MRCA Clade for *taxa* in a BioPython Phylo tree.

    Taxa not found in the tree are silently dropped.  If fewer than two
    valid taxa remain a warning is printed and ``None`` is returned.
    """
    tip_names = {tip.name for tip in tree.get_terminals()}
    valid = [t for t in taxa if t in tip_names]
    if len(valid) < 2:
        print(
            f"Warning: fewer than 2 valid taxa for MRCA lookup "
            f"(requested: {taxa}, valid: {valid}); skipping.",
            file=sys.stderr,
        )
        return None
    return tree.common_ancestor(valid)


def get_clade_tip_ids(clade) -> set:
    """Return ``{id(tip) for tip in clade.get_terminals()}``."""
    return {id(tip) for tip in clade.get_terminals()}


def get_clade_branch_ids(tree, clade, parent_map) -> set:
    """Return the set of ``id(node)`` for the MRCA and all its descendants."""
    ids = set()
    for node in clade.find_clades(order="preorder"):
        ids.add(id(node))
    return ids


# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------


def draw_range_rect(ax, tree, clade, color, node_x, node_y, alpha=0.15):
    """Draw a coloured rectangle behind *clade* in rectangular mode."""
    from matplotlib.patches import Rectangle

    tips = list(clade.get_terminals())
    if not tips:
        return

    tip_xs = [node_x[id(t)] for t in tips if id(t) in node_x]
    tip_ys = [node_y[id(t)] for t in tips if id(t) in node_y]
    if not tip_xs or not tip_ys:
        return

    mrca_x = node_x.get(id(clade), min(tip_xs))

    x_min = mrca_x
    x_max = max(tip_xs)
    # Small padding so the rectangle extends slightly past the tips
    x_pad = (x_max - x_min) * 0.05 if x_max > x_min else 0.02
    x_max += x_pad

    y_min = min(tip_ys) - 0.4
    y_max = max(tip_ys) + 0.4

    width = x_max - x_min
    height = y_max - y_min

    ax.add_patch(
        Rectangle(
            (x_min, y_min),
            width,
            height,
            facecolor=color,
            alpha=alpha,
            edgecolor="none",
            zorder=0,
        )
    )


def draw_range_wedge(ax, tree, clade, color, coords, alpha=0.15):
    """Draw a coloured wedge behind *clade* in circular mode.

    *coords* should be a dict mapping ``id(node)`` to ``(x, y)`` pairs in
    Cartesian coordinates (the circular layout).
    """
    from matplotlib.patches import Wedge

    tips = list(clade.get_terminals())
    if not tips:
        return

    tip_coords = [
        coords[id(t)] for t in tips if id(t) in coords
    ]
    if len(tip_coords) < 2:
        return

    # Compute angles and radii from Cartesian coords
    angles = [atan2(y, x) for x, y in tip_coords]
    radii = [(x ** 2 + y ** 2) ** 0.5 for x, y in tip_coords]

    # Sort angles for range computation
    sorted_angles = sorted(angles)

    # Compute angular gap between consecutive sorted tips to determine padding
    if len(sorted_angles) >= 2:
        gaps = [
            sorted_angles[i + 1] - sorted_angles[i]
            for i in range(len(sorted_angles) - 1)
        ]
        min_gap = min(gaps) if gaps else 0.05
        pad = min_gap * 0.5
    else:
        pad = 0.05

    # Check if the clade wraps around (largest gap > pi)
    # Compute the complement gap (from last to first going through 2*pi)
    complement_gap = (2 * pi) - (sorted_angles[-1] - sorted_angles[0])
    all_gaps = [
        sorted_angles[i + 1] - sorted_angles[i]
        for i in range(len(sorted_angles) - 1)
    ]
    all_gaps.append(complement_gap)

    # The largest gap is where the clade does NOT span
    max_gap_idx = all_gaps.index(max(all_gaps))
    if max_gap_idx < len(sorted_angles) - 1:
        angle_min = sorted_angles[max_gap_idx + 1] - pad
        angle_max = sorted_angles[max_gap_idx] + pad
        # This means we need to wrap
        if angle_max < angle_min:
            angle_max += 2 * pi
    else:
        # The biggest gap is the complement gap -> simple contiguous range
        angle_min = sorted_angles[0] - pad
        angle_max = sorted_angles[-1] + pad

    # Radii
    mrca_coord = coords.get(id(clade))
    if mrca_coord is not None:
        r_inner = (mrca_coord[0] ** 2 + mrca_coord[1] ** 2) ** 0.5
    else:
        r_inner = min(radii)

    r_outer = max(radii)
    r_pad = (r_outer - r_inner) * 0.05 if r_outer > r_inner else 0.02
    r_outer += r_pad

    ax.add_patch(
        Wedge(
            (0, 0),
            r_outer,
            degrees(angle_min),
            degrees(angle_max),
            width=r_outer - r_inner,
            facecolor=color,
            alpha=alpha,
            edgecolor="none",
            zorder=0,
        )
    )
