"""Shared circular (radial/fan) phylogram layout utilities.

Functions for computing polar coordinates from a phylogenetic tree and
drawing branches, arcs, tip labels, gradients, and multi-segment branches
on a circular layout.
"""

import math


_ARC_FRACTIONS = tuple(idx / 60 for idx in range(61))
_ARC_FRACTIONS_ARRAY = None


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()


def _arc_fractions_array():
    global _ARC_FRACTIONS_ARRAY
    if _ARC_FRACTIONS_ARRAY is None:
        _ARC_FRACTIONS_ARRAY = np.asarray(_ARC_FRACTIONS, dtype=float)
    return _ARC_FRACTIONS_ARRAY


# ---------------------------------------------------------------------------
# Coordinate computation
# ---------------------------------------------------------------------------

def compute_circular_coords(
    tree,
    node_x,
    parent_map,
    preorder_clades=None,
    terminal_clades=None,
):
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
    if preorder_clades is not None and terminal_clades is not None:
        coords = _compute_circular_coords_from_clades(
            tree,
            node_x,
            preorder_clades,
            terminal_clades,
        )
        if coords is not None:
            return coords

    try:
        root = tree.root
        root.clades
    except AttributeError:
        return _compute_circular_coords_legacy(tree, node_x, parent_map)

    clades = []
    stack = [root]
    while stack:
        clade = stack.pop()
        children = getattr(clade, "clades", None)
        if not isinstance(children, list):
            return _compute_circular_coords_legacy(tree, node_x, parent_map)
        clades.append(clade)
        if children:
            stack.extend(children)

    node_tip_range: dict[int, tuple[int, int]] = {}
    n_tips = 0

    for clade in reversed(clades):
        cid = id(clade)
        children = clade.clades
        if children:
            child_min, child_max = node_tip_range[id(children[0])]
            for child in children[1:]:
                cmin, cmax = node_tip_range[id(child)]
                if cmin < child_min:
                    child_min = cmin
                if cmax > child_max:
                    child_max = cmax
            node_tip_range[cid] = (child_min, child_max)
        else:
            tip_index = n_tips
            n_tips += 1
            node_tip_range[cid] = (tip_index, tip_index)

    coords: dict[int, dict[str, float]] = {}
    stack = [root]
    while stack:
        clade = stack.pop()
        cid = id(clade)
        radius = node_x[cid]
        if clade == root:
            angle = 0.0
        else:
            mn, mx = node_tip_range[cid]
            angle = math.pi * (mn + mx) / n_tips
        coords[cid] = {
            "x": radius * math.cos(angle),
            "y": radius * math.sin(angle),
            "angle": angle,
            "radius": radius,
        }

        children = clade.clades
        if children:
            stack.extend(reversed(children))

    return coords


def _compute_circular_coords_from_clades(
    tree,
    node_x,
    preorder_clades,
    terminal_clades,
):
    try:
        root = tree.root
        root.clades
    except AttributeError:
        return None

    n_tips = len(terminal_clades)
    if n_tips == 0:
        return None

    node_tip_range: dict[int, tuple[int, int]] = {
        id(tip): (idx, idx) for idx, tip in enumerate(terminal_clades)
    }

    try:
        for clade in reversed(preorder_clades):
            children = clade.clades
            cid = id(clade)
            if not children:
                if cid not in node_tip_range:
                    return None
                continue

            child_min, child_max = node_tip_range[id(children[0])]
            for child in children[1:]:
                cmin, cmax = node_tip_range[id(child)]
                if cmin < child_min:
                    child_min = cmin
                if cmax > child_max:
                    child_max = cmax
            node_tip_range[cid] = (child_min, child_max)
    except (AttributeError, KeyError, TypeError):
        return None

    coords: dict[int, dict[str, float]] = {}
    try:
        for clade in preorder_clades:
            cid = id(clade)
            radius = node_x[cid]
            if clade is root:
                angle = 0.0
            else:
                mn, mx = node_tip_range[cid]
                angle = math.pi * (mn + mx) / n_tips
            coords[cid] = {
                "x": radius * math.cos(angle),
                "y": radius * math.sin(angle),
                "angle": angle,
                "radius": radius,
            }
    except KeyError:
        return None

    return coords


def _compute_circular_coords_legacy(tree, node_x, parent_map):
    tips = tree.get_terminals()
    n_tips = len(tips)

    # Assign each tip an evenly-spaced angle.
    tip_angles: dict[int, float] = {}
    for i, tip in enumerate(tips):
        tip_angles[id(tip)] = 2.0 * math.pi * i / n_tips

    # We also need to know which tip indices belong to each internal node
    # so that we can compute the angular *range* of its descendants.
    # Walk postorder; for each node store (min_index, max_index) among tips.
    tip_index = {id(t): i for i, t in enumerate(tips)}
    node_tip_range: dict[int, tuple[int, int]] = {}

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
    node_angle: dict[int, float] = {}
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
    coords: dict[int, dict[str, float]] = {}
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


def _terminal_clades(tree):
    direct_terminals = _terminal_clades_direct(tree)
    if direct_terminals is not None:
        return direct_terminals
    return list(tree.get_terminals())


def _terminal_clades_direct(tree):
    try:
        root = tree.root
        root.clades
    except AttributeError:
        return None

    terminals = []
    stack = [root]
    pop = stack.pop
    append = stack.append
    append_terminal = terminals.append
    try:
        while stack:
            clade = pop()
            children = clade.clades
            if children:
                child_count = len(children)
                if child_count == 2:
                    append(children[1])
                    append(children[0])
                else:
                    for index in range(child_count - 1, -1, -1):
                        append(children[index])
            else:
                append_terminal(clade)
    except AttributeError:
        return None
    return terminals


def _preorder_clades(tree):
    direct_clades = _preorder_clades_direct(tree)
    if direct_clades is not None:
        return direct_clades
    return list(tree.find_clades(order="preorder"))


def _preorder_clades_direct(tree):
    try:
        root = tree.root
        root.clades
    except AttributeError:
        return None

    clades = []
    stack = [root]
    pop = stack.pop
    append = stack.append
    append_clade = clades.append
    try:
        while stack:
            clade = pop()
            append_clade(clade)
            children = clade.clades
            if children:
                child_count = len(children)
                if child_count == 2:
                    append(children[1])
                    append(children[0])
                else:
                    for index in range(child_count - 1, -1, -1):
                        append(children[index])
    except AttributeError:
        return None
    return clades


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

    angles = start + diff * _arc_fractions_array()
    xs = cx + radius * np.cos(angles)
    ys = cy + radius * np.sin(angles)
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
    if hasattr(ax, "add_collection"):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            pass
        else:
            _draw_circular_branches_collections(
                ax,
                root,
                coords,
                parent_map,
                color,
                lw,
            )
            return

    plot = ax.plot
    tau = 2.0 * math.pi
    arc_fractions = _arc_fractions_array()

    try:
        root = tree.root
        root.clades
    except AttributeError:
        clades = _preorder_clades(tree)
    else:
        clades = None

    if clades is None:
        stack = [root]
        pop = stack.pop
        append = stack.append
        while stack:
            clade = pop()
            children = clade.clades
            if children:
                child_count = len(children)
                if child_count == 2:
                    append(children[1])
                    append(children[0])
                else:
                    for index in range(child_count - 1, -1, -1):
                        append(children[index])

            cid = id(clade)
            if clade is not root:
                parent = parent_map[cid]
                pc = coords[id(parent)]
                cc = coords[cid]
                # Radial segment at child's angle.
                r_parent = pc["radius"]
                r_child = cc["radius"]
                x1 = cc["x"]
                y1 = cc["y"]
                if r_child:
                    scale = r_parent / r_child
                    x0 = x1 * scale
                    y0 = y1 * scale
                else:
                    x0 = 0.0
                    y0 = 0.0
                plot([x0, x1], [y0, y1], color=color, linewidth=lw,
                     solid_capstyle="round")

            if len(children) < 2:
                continue
            pc = coords[cid]

            start_a = None
            end_a = None
            for child in children:
                angle = coords[id(child)]["angle"]
                if start_a is None:
                    start_a = angle
                    end_a = angle
                else:
                    if angle < start_a:
                        start_a = angle
                    if angle > end_a:
                        end_a = angle
            span = (end_a - start_a) % tau
            if span > math.pi:
                # Swap so arc takes the long way (the actual subtended arc).
                start_a, end_a = end_a, start_a

            start = start_a % tau
            end = end_a % tau
            diff = (end - start) % tau
            if diff > math.pi:
                diff = diff - tau
            angles = start + diff * arc_fractions
            radius = pc["radius"]
            xs = radius * np.cos(angles)
            ys = radius * np.sin(angles)
            plot(xs, ys, color=color, linewidth=lw, solid_capstyle="round")
        return

    for clade in clades:
        cid = id(clade)
        if clade is not root:
            parent = parent_map[cid]
            pc = coords[id(parent)]
            cc = coords[cid]
            # Radial segment at child's angle.
            r_parent = pc["radius"]
            r_child = cc["radius"]
            x1 = cc["x"]
            y1 = cc["y"]
            if r_child:
                scale = r_parent / r_child
                x0 = x1 * scale
                y0 = y1 * scale
            else:
                x0 = 0.0
                y0 = 0.0
            plot([x0, x1], [y0, y1], color=color, linewidth=lw,
                 solid_capstyle="round")

        children = clade.clades
        if not children:
            continue
        pc = coords[id(clade)]

        start_a = None
        end_a = None
        child_count = 0
        for child in children:
            angle = coords[id(child)]["angle"]
            child_count += 1
            if start_a is None:
                start_a = angle
                end_a = angle
            else:
                if angle < start_a:
                    start_a = angle
                if angle > end_a:
                    end_a = angle
        if child_count < 2:
            continue
        span = (end_a - start_a) % (2.0 * math.pi)
        if span > math.pi:
            # Swap so arc takes the long way (the actual subtended arc).
            start_a, end_a = end_a, start_a
        _draw_arc(ax, 0.0, 0.0, pc["radius"], start_a, end_a, color, lw)


def _draw_circular_branches_collections(
    ax,
    root,
    coords,
    parent_map,
    color,
    lw,
) -> None:
    from matplotlib.collections import LineCollection

    tau = 2.0 * math.pi
    arc_fractions = _arc_fractions_array()
    radial_segments = []
    arc_segments = []
    stack = [root]
    coords_get = coords.__getitem__
    parent_get = parent_map.__getitem__
    id_ = id

    while stack:
        clade = stack.pop()
        children = clade.clades
        if children:
            stack.extend(reversed(children))

        cid = id_(clade)
        if clade is not root:
            pc = coords_get(id_(parent_get(cid)))
            cc = coords_get(cid)
            r_parent = pc["radius"]
            r_child = cc["radius"]
            x1 = cc["x"]
            y1 = cc["y"]
            if r_child:
                scale = r_parent / r_child
                x0 = x1 * scale
                y0 = y1 * scale
            else:
                x0 = 0.0
                y0 = 0.0
            radial_segments.append(((x0, y0), (x1, y1)))

        if len(children) < 2:
            continue

        start_a = None
        end_a = None
        for child in children:
            angle = coords_get(id_(child))["angle"]
            if start_a is None:
                start_a = angle
                end_a = angle
            else:
                if angle < start_a:
                    start_a = angle
                if angle > end_a:
                    end_a = angle

        span = (end_a - start_a) % tau
        if span > math.pi:
            start_a, end_a = end_a, start_a

        start = start_a % tau
        end = end_a % tau
        diff = (end - start) % tau
        if diff > math.pi:
            diff -= tau
        angles = start + diff * arc_fractions
        radius = coords_get(cid)["radius"]
        arc_segments.append(
            np.column_stack((radius * np.cos(angles), radius * np.sin(angles)))
        )

    if radial_segments:
        ax.add_collection(
            LineCollection(
                radial_segments,
                colors=color,
                linewidths=lw,
                capstyle="round",
            ),
            autolim=True,
        )
    if arc_segments:
        ax.add_collection(
            LineCollection(
                arc_segments,
                colors=color,
                linewidths=lw,
                capstyle="round",
            ),
            autolim=True,
        )
    ax.autoscale_view()


def draw_circular_tip_labels(ax, tree, coords, fontsize=9, offset=0.03):
    """Place tip labels around the circular layout.

    Labels on the right half (-90 to 90 degrees) are left-aligned with
    rotation following the angle.  Labels on the left half are
    right-aligned and rotated by 180 degrees so text always reads
    left-to-right.
    """
    text = ax.text
    coords_get = coords.__getitem__
    id_ = id
    cos = math.cos
    sin = math.sin
    degrees = math.degrees

    for tip in _terminal_clades(tree):
        tc = coords_get(id_(tip))
        angle = tc["angle"]
        r = tc["radius"] + offset
        x = r * cos(angle)
        y = r * sin(angle)

        deg = degrees(angle)
        # Normalise to (-180, 180]
        deg = ((deg + 180) % 360) - 180

        if -90 <= deg <= 90:
            ha = "left"
            rotation = deg
        else:
            ha = "right"
            rotation = deg + 180

        label = tip.name if tip.name else ""
        text(x, y, label, fontsize=fontsize, ha=ha, va="center",
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


def draw_circular_colored_arcs(ax, arcs, lw=1.5):
    """Draw many coloured circular arcs as one collection when possible.

    ``arcs`` contains ``(cx, cy, radius, start_angle, end_angle, color)`` tuples.
    """
    if not arcs:
        return

    if not hasattr(ax, "add_collection"):
        for cx, cy, radius, start_angle, end_angle, color in arcs:
            draw_circular_colored_arc(
                ax, cx, cy, radius, start_angle, end_angle, color, lw=lw,
            )
        return

    from matplotlib.collections import LineCollection

    arc_fractions = _arc_fractions_array()
    arc_segments = []
    colors = []
    for cx, cy, radius, start_angle, end_angle, color in arcs:
        start = start_angle % (2.0 * math.pi)
        end = end_angle % (2.0 * math.pi)
        diff = (end - start) % (2.0 * math.pi)
        if diff > math.pi:
            diff = diff - 2.0 * math.pi

        angles = start + diff * arc_fractions
        xs = cx + radius * np.cos(angles)
        ys = cy + radius * np.sin(angles)
        arc_segments.append(np.column_stack((xs, ys)))
        colors.append(color)

    ax.add_collection(
        LineCollection(
            arc_segments,
            colors=colors,
            linewidths=lw,
            capstyle="round",
        ),
        autolim=True,
    )
    ax.autoscale_view()


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
    cos_angle = math.cos(angle)
    sin_angle = math.sin(angle)
    denominator = max(n_points - 1, 1)
    radius_delta = r_c - r_p
    points: list[tuple[float, float, float]] = []
    for i in range(n_points):
        r = r_p + radius_delta * (i / denominator)
        points.append((r * cos_angle, r * sin_angle, angle))
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
    cos_angle = math.cos(angle)
    sin_angle = math.sin(angle)

    if hasattr(ax, "add_collection"):
        from matplotlib.collections import LineCollection

        segments = []
        colors = []
        radius_delta = r_c - r_p
        for i in range(n_seg):
            t0 = i / n_seg
            t1 = (i + 1) / n_seg
            r0 = r_p + radius_delta * t0
            r1 = r_p + radius_delta * t1
            val = parent_val + (child_val - parent_val) * (t0 + t1) / 2.0
            normed = (val - vmin) / norm_range
            normed = max(0.0, min(1.0, normed))
            segments.append((
                (r0 * cos_angle, r0 * sin_angle),
                (r1 * cos_angle, r1 * sin_angle),
            ))
            colors.append(cmap(normed))
        ax.add_collection(
            LineCollection(
                segments,
                colors=colors,
                linewidths=lw,
                capstyle="butt",
            ),
            autolim=True,
        )
        ax.autoscale_view()
        return

    for i in range(n_seg):
        t0 = i / n_seg
        t1 = (i + 1) / n_seg
        r0 = r_p + (r_c - r_p) * t0
        r1 = r_p + (r_c - r_p) * t1
        val = parent_val + (child_val - parent_val) * (t0 + t1) / 2.0
        normed = (val - vmin) / norm_range
        normed = max(0.0, min(1.0, normed))
        seg_color = cmap(normed)
        x0 = r0 * cos_angle
        y0 = r0 * sin_angle
        x1 = r1 * cos_angle
        y1 = r1 * sin_angle
        ax.plot([x0, x1], [y0, y1], color=seg_color, linewidth=lw,
                solid_capstyle="butt")


def draw_circular_gradient_branches(ax, branch_data, cmap, vmin, vmax, lw=1.5):
    """Draw many radial gradient branches as one collection when possible.

    ``branch_data`` contains ``(parent_coords, child_coords, parent_val,
    child_val)`` tuples matching :func:`draw_circular_gradient_branch`.
    """
    if not branch_data:
        return

    if not hasattr(ax, "add_collection"):
        for parent_coords, child_coords, parent_val, child_val in branch_data:
            draw_circular_gradient_branch(
                ax, parent_coords, child_coords,
                cmap, vmin, vmax, parent_val, child_val, lw=lw,
            )
        return

    from matplotlib.collections import LineCollection

    n_seg = 30
    norm_range = vmax - vmin if vmax != vmin else 1.0
    segments = []
    colors = []

    for parent_coords, child_coords, parent_val, child_val in branch_data:
        angle = child_coords["angle"]
        r_p = parent_coords["radius"]
        r_c = child_coords["radius"]
        radius_delta = r_c - r_p
        cos_angle = math.cos(angle)
        sin_angle = math.sin(angle)

        for i in range(n_seg):
            t0 = i / n_seg
            t1 = (i + 1) / n_seg
            r0 = r_p + radius_delta * t0
            r1 = r_p + radius_delta * t1
            val = parent_val + (child_val - parent_val) * (t0 + t1) / 2.0
            normed = (val - vmin) / norm_range
            normed = max(0.0, min(1.0, normed))
            segments.append((
                (r0 * cos_angle, r0 * sin_angle),
                (r1 * cos_angle, r1 * sin_angle),
            ))
            colors.append(cmap(normed))

    ax.add_collection(
        LineCollection(
            segments,
            colors=colors,
            linewidths=lw,
            capstyle="butt",
        ),
        autolim=True,
    )
    ax.autoscale_view()


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
    cos_angle = math.cos(angle)
    sin_angle = math.sin(angle)

    if hasattr(ax, "add_collection"):
        from matplotlib.collections import LineCollection

        branch_segments = []
        colors = []
        radius_delta = r_c - r_p
        for start_frac, end_frac, state in segments:
            r0 = r_p + radius_delta * start_frac
            r1 = r_p + radius_delta * end_frac
            branch_segments.append((
                (r0 * cos_angle, r0 * sin_angle),
                (r1 * cos_angle, r1 * sin_angle),
            ))
            colors.append(state_colors[state])
        if branch_segments:
            ax.add_collection(
                LineCollection(
                    branch_segments,
                    colors=colors,
                    linewidths=lw,
                    capstyle="butt",
                ),
                autolim=True,
            )
            ax.autoscale_view()
        return

    for start_frac, end_frac, state in segments:
        r0 = r_p + (r_c - r_p) * start_frac
        r1 = r_p + (r_c - r_p) * end_frac
        x0 = r0 * cos_angle
        y0 = r0 * sin_angle
        x1 = r1 * cos_angle
        y1 = r1 * sin_angle
        ax.plot([x0, x1], [y0, y1], color=state_colors[state],
                linewidth=lw, solid_capstyle="butt")
