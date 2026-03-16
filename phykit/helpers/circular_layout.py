"""Shared circular (radial/fan) phylogram layout utilities.

Functions for computing polar coordinates from a phylogenetic tree and
drawing branches, arcs, tip labels, gradients, and multi-segment branches
on a circular layout.
"""

import math
from typing import Any, Dict, List, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Coordinate computation
# ---------------------------------------------------------------------------

def compute_circular_coords(tree, node_x, parent_map):
    """Compute Cartesian coordinates for a circular phylogram layout.

    Each tip is evenly spaced around the circle.  Internal node angles are
    the midpoint of the angular range of their descendant tips (NOT the
    circular mean, which fails when descendants span the full circle).
    The root gets angle = 0 (irrelevant since its radius is 0).

    Parameters
    ----------
    tree : Bio.Phylo tree
        The phylogenetic tree.
    node_x : dict
        Mapping ``id(clade) -> float`` giving the radial distance (reuses
        existing phylogram / cladogram x-values).
    parent_map : dict
        Mapping ``id(child) -> parent_clade``.

    Returns
    -------
    dict
        ``id(clade) -> {"x": float, "y": float, "angle": float, "radius": float}``
    """
    tips = tree.get_terminals()
    n_tips = len(tips)

    # Assign each tip an evenly-spaced angle.
    tip_angles: Dict[int, float] = {}
    for i, tip in enumerate(tips):
        tip_angles[id(tip)] = 2.0 * math.pi * i / n_tips

    # We also need to know which tip indices belong to each internal node
    # so that we can compute the angular *range* of its descendants.
    # Walk postorder; for each node store (min_index, max_index) among tips.
    tip_index = {id(t): i for i, t in enumerate(tips)}
    node_tip_range: Dict[int, Tuple[int, int]] = {}

    for clade in tree.find_clades(order="postorder"):
        cid = id(clade)
        if clade.is_terminal():
            idx = tip_index[cid]
            node_tip_range[cid] = (idx, idx)
        else:
            child_mins = []
            child_maxs = []
            for child in clade.clades:
                cmin, cmax = node_tip_range[id(child)]
                child_mins.append(cmin)
                child_maxs.append(cmax)
            node_tip_range[cid] = (min(child_mins), max(child_maxs))

    # Build angle dict for every node.
    node_angle: Dict[int, float] = {}
    root = tree.root
    for clade in tree.find_clades(order="postorder"):
        cid = id(clade)
        if clade.is_terminal():
            node_angle[cid] = tip_angles[cid]
        elif clade == root:
            node_angle[cid] = 0.0
        else:
            mn, mx = node_tip_range[cid]
            a_min = 2.0 * math.pi * mn / n_tips
            a_max = 2.0 * math.pi * mx / n_tips
            node_angle[cid] = (a_min + a_max) / 2.0

    # Convert to Cartesian.
    coords: Dict[int, Dict[str, float]] = {}
    for clade in tree.find_clades(order="preorder"):
        cid = id(clade)
        radius = node_x[cid]
        angle = node_angle[cid]
        coords[cid] = {
            "x": radius * math.cos(angle),
            "y": radius * math.sin(angle),
            "angle": angle,
            "radius": radius,
        }
    return coords


# ---------------------------------------------------------------------------
# Arc helper (internal)
# ---------------------------------------------------------------------------

def _draw_arc(ax, cx, cy, radius, start_angle, end_angle, color, lw):
    """Draw a circular arc as a polyline.

    Sweeps from *start_angle* to *end_angle* at the given *radius* around
    the centre (*cx*, *cy*).  Handles wrap-around so the arc always takes
    the shorter path unless the caller has already adjusted the angles.
    """
    # Normalise into [0, 2pi)
    start = start_angle % (2.0 * math.pi)
    end = end_angle % (2.0 * math.pi)

    # Choose the shorter sweep direction.
    diff = (end - start) % (2.0 * math.pi)
    if diff > math.pi:
        # Sweep the other way.
        diff = diff - 2.0 * math.pi

    n_pts = 61
    angles = [start + diff * t / (n_pts - 1) for t in range(n_pts)]
    xs = [cx + radius * math.cos(a) for a in angles]
    ys = [cy + radius * math.sin(a) for a in angles]
    ax.plot(xs, ys, color=color, linewidth=lw, solid_capstyle="round")


# ---------------------------------------------------------------------------
# Drawing functions
# ---------------------------------------------------------------------------

def draw_circular_branches(ax, tree, coords, parent_map,
                           color="#333", lw=1.5):
    """Draw all branches of a circular phylogram.

    For each non-root clade a *radial* line segment is drawn from the
    parent radius to the child radius at the child's angle.  For each
    internal node an *arc* is drawn at the parent's radius connecting the
    children's angles.
    """
    root = tree.root

    for clade in tree.find_clades(order="preorder"):
        cid = id(clade)
        if clade == root:
            continue
        parent = parent_map[cid]
        pc = coords[id(parent)]
        cc = coords[cid]
        # Radial segment at child's angle.
        r_parent = pc["radius"]
        r_child = cc["radius"]
        angle = cc["angle"]
        x0 = r_parent * math.cos(angle)
        y0 = r_parent * math.sin(angle)
        x1 = r_child * math.cos(angle)
        y1 = r_child * math.sin(angle)
        ax.plot([x0, x1], [y0, y1], color=color, linewidth=lw,
                solid_capstyle="round")

    # Arcs for internal nodes.
    for clade in tree.find_clades(order="preorder"):
        if clade.is_terminal() or not clade.clades:
            continue
        pc = coords[id(clade)]
        child_angles = [coords[id(ch)]["angle"] for ch in clade.clades]
        if len(child_angles) < 2:
            continue
        start_a = min(child_angles)
        end_a = max(child_angles)
        span = (end_a - start_a) % (2.0 * math.pi)
        if span > math.pi:
            # Swap so arc takes the long way (the actual subtended arc).
            start_a, end_a = end_a, start_a
        _draw_arc(ax, 0.0, 0.0, pc["radius"], start_a, end_a, color, lw)


def draw_circular_tip_labels(ax, tree, coords, fontsize=9, offset=0.03):
    """Place tip labels around the circular layout.

    Labels on the right half (-90 to 90 degrees) are left-aligned with
    rotation following the angle.  Labels on the left half are
    right-aligned and rotated by 180 degrees so text always reads
    left-to-right.
    """
    for tip in tree.get_terminals():
        tc = coords[id(tip)]
        angle = tc["angle"]
        r = tc["radius"] + offset
        x = r * math.cos(angle)
        y = r * math.sin(angle)

        deg = math.degrees(angle)
        # Normalise to (-180, 180]
        deg = ((deg + 180) % 360) - 180

        if -90 <= deg <= 90:
            ha = "left"
            rotation = deg
        else:
            ha = "right"
            rotation = deg + 180

        label = tip.name if tip.name else ""
        ax.text(x, y, label, fontsize=fontsize, ha=ha, va="center",
                rotation=rotation, rotation_mode="anchor")


def draw_circular_colored_branch(ax, parent_coords, child_coords,
                                 color, lw=1.5):
    """Draw a single radial segment in a specific colour.

    Goes from ``parent_radius`` at the child's angle to ``child_radius``
    at the child's angle.
    """
    angle = child_coords["angle"]
    r_p = parent_coords["radius"]
    r_c = child_coords["radius"]
    x0 = r_p * math.cos(angle)
    y0 = r_p * math.sin(angle)
    x1 = r_c * math.cos(angle)
    y1 = r_c * math.sin(angle)
    ax.plot([x0, x1], [y0, y1], color=color, linewidth=lw,
            solid_capstyle="round")


def draw_circular_colored_arc(ax, cx, cy, radius, start_angle, end_angle,
                              color, lw=1.5):
    """Public wrapper around ``_draw_arc`` for coloured arcs."""
    _draw_arc(ax, cx, cy, radius, start_angle, end_angle, color, lw)


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

def circular_branch_points(parent_coords, child_coords, n_points):
    """Return *n_points* evenly spaced ``(x, y, angle)`` tuples along a
    radial segment from parent to child.
    """
    angle = child_coords["angle"]
    r_p = parent_coords["radius"]
    r_c = child_coords["radius"]
    points: List[Tuple[float, float, float]] = []
    for i in range(n_points):
        t = i / max(n_points - 1, 1)
        r = r_p + (r_c - r_p) * t
        points.append((r * math.cos(angle), r * math.sin(angle), angle))
    return points


def radial_offset(angle, distance):
    """Return ``(dx, dy)`` for a displacement of *distance* at *angle*."""
    return (math.cos(angle) * distance, math.sin(angle) * distance)


# ---------------------------------------------------------------------------
# Gradient and multi-segment branches
# ---------------------------------------------------------------------------

def draw_circular_gradient_branch(ax, parent_coords, child_coords,
                                  cmap, vmin, vmax, parent_val, child_val,
                                  lw=1.5):
    """Draw a radial segment coloured by a continuous gradient.

    The segment is divided into ~30 sub-segments; each is coloured by
    linearly interpolating *parent_val* to *child_val* through *cmap*.
    """
    n_seg = 30
    angle = child_coords["angle"]
    r_p = parent_coords["radius"]
    r_c = child_coords["radius"]
    norm_range = vmax - vmin if vmax != vmin else 1.0

    for i in range(n_seg):
        t0 = i / n_seg
        t1 = (i + 1) / n_seg
        r0 = r_p + (r_c - r_p) * t0
        r1 = r_p + (r_c - r_p) * t1
        val = parent_val + (child_val - parent_val) * (t0 + t1) / 2.0
        normed = (val - vmin) / norm_range
        normed = max(0.0, min(1.0, normed))
        seg_color = cmap(normed)
        x0 = r0 * math.cos(angle)
        y0 = r0 * math.sin(angle)
        x1 = r1 * math.cos(angle)
        y1 = r1 * math.sin(angle)
        ax.plot([x0, x1], [y0, y1], color=seg_color, linewidth=lw,
                solid_capstyle="butt")


def draw_circular_multi_segment_branch(ax, parent_coords, child_coords,
                                       segments, state_colors, lw=1.5):
    """Draw a radial segment split into coloured sub-segments.

    Parameters
    ----------
    segments : list of (start_frac, end_frac, state)
        Each tuple gives the fractional start/end along the branch and
        the discrete state label.
    state_colors : dict
        Mapping ``state -> colour``.
    """
    angle = child_coords["angle"]
    r_p = parent_coords["radius"]
    r_c = child_coords["radius"]

    for start_frac, end_frac, state in segments:
        r0 = r_p + (r_c - r_p) * start_frac
        r1 = r_p + (r_c - r_p) * end_frac
        x0 = r0 * math.cos(angle)
        y0 = r0 * math.sin(angle)
        x1 = r1 * math.cos(angle)
        y1 = r1 * math.sin(angle)
        ax.plot([x0, x1], [y0, y1], color=state_colors[state],
                linewidth=lw, solid_capstyle="butt")
