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
    labels: dict = {}
    ranges: list = []
    clades: list = []
    ranges_append = ranges.append
    clades_append = clades.append

    try:
        with open(path) as fh:
            for line in fh:
                stripped = line.strip()
                if not stripped or stripped[0] == "#":
                    continue

                parts = stripped.split("\t", 4)
                if len(parts) < 3:
                    continue

                field2 = parts[1].strip()
                if field2 not in ("label", "range", "clade"):
                    field2 = field2.lower()

                if field2 == "label":
                    labels[parts[0].strip()] = parts[2].strip()
                    continue
                if field2 == "range" or field2 == "clade":
                    field1 = parts[0].strip()
                    field3 = parts[2].strip()
                    field4 = parts[3].strip() if len(parts) > 3 else None
                    taxa = [t.strip() for t in field1.split(",")]
                    if field2 == "range":
                        ranges_append((taxa, field3, field4))
                    else:
                        clades_append((taxa, field3, field4))
    except FileNotFoundError:
        raise PhykitUserError([f"Color file not found: {path}"])

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
    valid = _valid_mrca_taxa(getattr(tree, "root", tree), taxa)
    if len(valid) < 2:
        print(
            f"Warning: fewer than 2 valid taxa for MRCA lookup "
            f"(requested: {taxa}, valid: {valid}); skipping.",
            file=sys.stderr,
        )
        return None
    return tree.common_ancestor(valid)


def _valid_mrca_taxa(root, taxa):
    requested = set(taxa)
    if not requested:
        return []

    found = set()
    stack = [root]
    pop = stack.pop
    append = stack.append
    while stack:
        node = pop()
        try:
            children = node.clades
        except AttributeError:
            tip_names = {tip.name for tip in _terminal_clades(root)}
            return [t for t in taxa if t in tip_names]
        if not isinstance(children, list):
            tip_names = {tip.name for tip in _terminal_clades(root)}
            return [t for t in taxa if t in tip_names]
        if children:
            child_count = len(children)
            if child_count == 2:
                append(children[1])
                append(children[0])
            else:
                for index in range(child_count - 1, -1, -1):
                    append(children[index])
            continue

        name = getattr(node, "name", None)
        if name in requested:
            found.add(name)
            if len(found) == len(requested):
                break

    return [t for t in taxa if t in found]


def _terminal_clades(clade):
    try:
        clade.clades
    except AttributeError:
        return list(clade.get_terminals())

    terminals = []
    stack = [clade]
    pop = stack.pop
    append = stack.append
    append_terminal = terminals.append
    try:
        while stack:
            node = pop()
            children = node.clades
            if children:
                child_count = len(children)
                if child_count == 2:
                    append(children[1])
                    append(children[0])
                else:
                    for index in range(child_count - 1, -1, -1):
                        append(children[index])
            else:
                append_terminal(node)
    except AttributeError:
        return list(clade.get_terminals())
    return terminals


def get_clade_tip_ids(clade) -> set:
    """Return ``{id(tip) for tip in clade.get_terminals()}``."""
    try:
        clade.clades
    except AttributeError:
        return {id(tip) for tip in clade.get_terminals()}

    ids = set()
    stack = [clade]
    pop = stack.pop
    extend = stack.extend
    add = ids.add
    while stack:
        node = pop()
        try:
            children = node.clades
        except AttributeError:
            return {id(tip) for tip in clade.get_terminals()}
        if not isinstance(children, list):
            return {id(tip) for tip in clade.get_terminals()}
        if children:
            extend(children)
        else:
            add(id(node))
    return ids


def get_clade_branch_ids(tree, clade, parent_map) -> set:
    """Return the set of ``id(node)`` for the MRCA and all its descendants."""
    ids = set()
    stack = [clade]
    pop = stack.pop
    extend = stack.extend
    add = ids.add
    try:
        while stack:
            node = pop()
            add(id(node))
            children = node.clades
            if children:
                extend(children)
    except AttributeError:
        ids.clear()
        for node in clade.find_clades(order="preorder"):
            ids.add(id(node))
    return ids


def build_color_legend_handles(color_data):
    """Build matplotlib legend handles for labeled ranges and clades.

    Returns a list of Patch objects for any range or clade entry that
    has a label (field 4). Import matplotlib lazily.
    """
    from matplotlib.patches import Patch

    handles = []
    for taxa_list, color, label in color_data.get("ranges", []):
        if label:
            handles.append(
                Patch(facecolor=color, alpha=0.3, edgecolor="none", label=label)
            )
    for taxa_list, color, label in color_data.get("clades", []):
        if label:
            handles.append(
                Patch(facecolor=color, edgecolor="black", linewidth=0.5, label=label)
            )
    return handles


def apply_label_colors(ax, label_colors):
    """Apply parsed label colors to matching Matplotlib text artists."""
    if not label_colors:
        return 0

    pending_labels = set(label_colors)
    colored_count = 0
    for text_obj in ax.texts:
        text = text_obj.get_text()
        if text not in pending_labels:
            continue

        text_obj.set_color(label_colors[text])
        pending_labels.remove(text)
        colored_count += 1
        if not pending_labels:
            break
    return colored_count


def _range_wedge_angle_bounds(sorted_angles):
    two_pi = 2 * pi
    if len(sorted_angles) >= 2:
        min_gap = sorted_angles[1] - sorted_angles[0]
        max_gap = min_gap
        max_gap_idx = 0
        for idx in range(1, len(sorted_angles) - 1):
            gap = sorted_angles[idx + 1] - sorted_angles[idx]
            if gap < min_gap:
                min_gap = gap
            if gap > max_gap:
                max_gap = gap
                max_gap_idx = idx
        pad = min_gap * 0.5
        complement_gap = two_pi - (sorted_angles[-1] - sorted_angles[0])
        if complement_gap > max_gap:
            max_gap_idx = len(sorted_angles) - 1
    else:
        pad = 0.05
        max_gap_idx = len(sorted_angles) - 1

    if max_gap_idx < len(sorted_angles) - 1:
        angle_min = sorted_angles[max_gap_idx + 1] - pad
        angle_max = sorted_angles[max_gap_idx] + pad
        if angle_max < angle_min:
            angle_max += two_pi
    else:
        angle_min = sorted_angles[0] - pad
        angle_max = sorted_angles[-1] + pad
    return angle_min, angle_max


def _range_rect_tip_bounds(tips, node_x, node_y):
    have_x = False
    have_y = False
    min_x = max_x = 0.0
    min_y = max_y = 0.0

    for tip in tips:
        tip_id = id(tip)
        if tip_id in node_x:
            x = node_x[tip_id]
            if have_x:
                if x < min_x:
                    min_x = x
                elif x > max_x:
                    max_x = x
            else:
                min_x = max_x = x
                have_x = True

        if tip_id in node_y:
            y = node_y[tip_id]
            if have_y:
                if y < min_y:
                    min_y = y
                elif y > max_y:
                    max_y = y
            else:
                min_y = max_y = y
                have_y = True

    if not have_x or not have_y:
        return None
    return min_x, max_x, min_y, max_y


# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------


def draw_range_rect(ax, tree, clade, color, node_x, node_y, alpha=0.15):
    """Draw a coloured rectangle behind *clade* in rectangular mode."""
    from matplotlib.patches import Rectangle

    tips = _terminal_clades(clade)
    if not tips:
        return

    tip_bounds = _range_rect_tip_bounds(tips, node_x, node_y)
    if tip_bounds is None:
        return
    min_tip_x, x_max, y_min, y_max = tip_bounds

    mrca_x = node_x.get(id(clade), min_tip_x)

    x_min = mrca_x
    # Small padding so the rectangle extends slightly past the tips
    x_pad = (x_max - x_min) * 0.05 if x_max > x_min else 0.02
    x_max += x_pad

    y_min -= 0.4
    y_max += 0.4

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

    tips = _terminal_clades(clade)
    if not tips:
        return

    angles = []
    append_angle = angles.append
    have_radius = False
    min_radius = max_radius = 0.0
    coords_get = coords.get
    for tip in tips:
        coord = coords_get(id(tip))
        if coord is None:
            continue
        append_angle(coord["angle"])
        radius = coord["radius"]
        if have_radius:
            if radius < min_radius:
                min_radius = radius
            elif radius > max_radius:
                max_radius = radius
        else:
            min_radius = max_radius = radius
            have_radius = True

    if len(angles) < 2:
        return

    # Sort angles for range computation
    sorted_angles = sorted(angles)

    angle_min, angle_max = _range_wedge_angle_bounds(sorted_angles)

    # Radii
    mrca_coord = coords.get(id(clade))
    if mrca_coord is not None:
        r_inner = mrca_coord["radius"]
    else:
        r_inner = min_radius

    r_outer = max_radius
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
