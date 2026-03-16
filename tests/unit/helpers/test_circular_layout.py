"""Tests for phykit.helpers.circular_layout."""

import math
from io import StringIO

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest
from Bio import Phylo

from phykit.helpers.circular_layout import (
    circular_branch_points,
    compute_circular_coords,
    draw_circular_branches,
    draw_circular_colored_arc,
    draw_circular_colored_branch,
    draw_circular_gradient_branch,
    draw_circular_multi_segment_branch,
    draw_circular_tip_labels,
    radial_offset,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_tree(newick):
    """Return a Bio.Phylo tree from a Newick string."""
    return Phylo.read(StringIO(newick), "newick")


def _build_parent_map(tree):
    parent_map = {}
    for clade in tree.find_clades(order="preorder"):
        for child in clade.clades:
            parent_map[id(child)] = clade
    return parent_map


def _build_node_x_phylogram(tree):
    """Simple cumulative branch-length x-coordinates."""
    node_x = {}
    for clade in tree.find_clades(order="preorder"):
        if clade == tree.root:
            node_x[id(clade)] = 0.0
        else:
            parent_map = _build_parent_map(tree)
            parent = parent_map[id(clade)]
            bl = clade.branch_length if clade.branch_length else 0.0
            node_x[id(clade)] = node_x[id(parent)] + bl
    return node_x


def _build_node_x_cladogram(tree):
    """All tips at radius 1.0, internals at topological depth fraction."""
    parent_map = _build_parent_map(tree)
    depth = {}
    for clade in tree.find_clades(order="preorder"):
        if clade == tree.root:
            depth[id(clade)] = 0
        else:
            depth[id(clade)] = depth[id(parent_map[id(clade)])] + 1
    max_d = max(depth.values()) if depth else 1
    node_x = {}
    for clade in tree.find_clades(order="preorder"):
        if clade.is_terminal():
            node_x[id(clade)] = 1.0
        else:
            node_x[id(clade)] = depth[id(clade)] / max(max_d, 1)
    return node_x


# Reusable simple tree: ((A:1,B:1):1,(C:1,D:1):1,(E:1,F:1):1);
NEWICK_6 = "((A:1,B:1):1,(C:1,D:1):1,(E:1,F:1):1);"


# ---------------------------------------------------------------------------
# Tests: compute_circular_coords
# ---------------------------------------------------------------------------

class TestTipAngles:
    def test_tip_angles_evenly_spaced(self):
        tree = _make_tree(NEWICK_6)
        parent_map = _build_parent_map(tree)
        node_x = _build_node_x_phylogram(tree)
        coords = compute_circular_coords(tree, node_x, parent_map)

        tips = tree.get_terminals()
        assert len(tips) == 6
        expected_step = 2.0 * math.pi / 6.0
        for i, tip in enumerate(tips):
            expected_angle = expected_step * i
            assert coords[id(tip)]["angle"] == pytest.approx(expected_angle, abs=1e-9)

    def test_tip_angles_monotonic(self):
        tree = _make_tree(NEWICK_6)
        parent_map = _build_parent_map(tree)
        node_x = _build_node_x_phylogram(tree)
        coords = compute_circular_coords(tree, node_x, parent_map)

        tips = tree.get_terminals()
        angles = [coords[id(t)]["angle"] for t in tips]
        for i in range(len(angles) - 1):
            assert angles[i] < angles[i + 1]


class TestInternalNodeAngles:
    def test_internal_node_angle_between_children(self):
        tree = _make_tree("((A:1,B:1):1,C:1);")
        parent_map = _build_parent_map(tree)
        node_x = _build_node_x_phylogram(tree)
        coords = compute_circular_coords(tree, node_x, parent_map)

        tips = tree.get_terminals()
        # Tips: A(index 0), B(index 1), C(index 2)
        # Internal (A,B) should have angle = midpoint of A and B angles.
        a_A = coords[id(tips[0])]["angle"]
        a_B = coords[id(tips[1])]["angle"]
        # Find the internal node that is the parent of A.
        internal = parent_map[id(tips[0])]
        a_int = coords[id(internal)]["angle"]
        mid = (a_A + a_B) / 2.0
        assert a_int == pytest.approx(mid, abs=1e-9)

    def test_internal_angle_large_subtree(self):
        # 5 tips under one clade, 1 tip outside -> subtree spans > 180 deg.
        tree = _make_tree("((A:1,B:1,C:1,D:1,E:1):1,F:1);")
        parent_map = _build_parent_map(tree)
        node_x = _build_node_x_phylogram(tree)
        coords = compute_circular_coords(tree, node_x, parent_map)

        tips = tree.get_terminals()
        # The big internal node parents tips 0..4.
        big_internal = parent_map[id(tips[0])]
        a = coords[id(big_internal)]["angle"]
        # Midpoint of indices 0..4 out of 6 tips.
        a_min = 2.0 * math.pi * 0 / 6
        a_max = 2.0 * math.pi * 4 / 6
        expected = (a_min + a_max) / 2.0
        assert a == pytest.approx(expected, abs=1e-9)


class TestRootAndRadius:
    def test_root_at_center(self):
        tree = _make_tree(NEWICK_6)
        parent_map = _build_parent_map(tree)
        node_x = _build_node_x_phylogram(tree)
        # Root has radius 0 -> x,y ~ 0.
        coords = compute_circular_coords(tree, node_x, parent_map)
        rc = coords[id(tree.root)]
        assert rc["x"] == pytest.approx(0.0, abs=1e-12)
        assert rc["y"] == pytest.approx(0.0, abs=1e-12)

    def test_coords_on_circle(self):
        tree = _make_tree(NEWICK_6)
        parent_map = _build_parent_map(tree)
        node_x = _build_node_x_phylogram(tree)
        coords = compute_circular_coords(tree, node_x, parent_map)
        for clade in tree.find_clades():
            c = coords[id(clade)]
            r = c["radius"]
            dist = math.sqrt(c["x"] ** 2 + c["y"] ** 2)
            assert dist == pytest.approx(r, abs=1e-9)

    def test_cladogram_tips_same_radius(self):
        tree = _make_tree(NEWICK_6)
        parent_map = _build_parent_map(tree)
        node_x = _build_node_x_cladogram(tree)
        coords = compute_circular_coords(tree, node_x, parent_map)
        tips = tree.get_terminals()
        radii = [coords[id(t)]["radius"] for t in tips]
        for r in radii:
            assert r == pytest.approx(radii[0], abs=1e-9)


# ---------------------------------------------------------------------------
# Tests: circular_branch_points
# ---------------------------------------------------------------------------

class TestBranchPoints:
    def test_branch_points_count(self):
        pc = {"angle": 0.5, "radius": 1.0}
        cc = {"angle": 0.5, "radius": 2.0}
        pts = circular_branch_points(pc, cc, 10)
        assert len(pts) == 10

    def test_branch_points_on_radial_line(self):
        angle = 1.2
        pc = {"angle": angle, "radius": 0.5}
        cc = {"angle": angle, "radius": 3.0}
        pts = circular_branch_points(pc, cc, 20)
        for x, y, a in pts:
            assert a == pytest.approx(angle, abs=1e-12)
            # Check that atan2 gives the same angle.
            r = math.sqrt(x ** 2 + y ** 2)
            if r > 1e-12:
                computed_angle = math.atan2(y, x)
                assert computed_angle == pytest.approx(angle, abs=1e-6)


# ---------------------------------------------------------------------------
# Tests: radial_offset
# ---------------------------------------------------------------------------

class TestRadialOffset:
    def test_radial_offset_direction(self):
        d = 2.5
        dx, dy = radial_offset(0.0, d)
        assert dx == pytest.approx(d, abs=1e-12)
        assert dy == pytest.approx(0.0, abs=1e-12)

        dx, dy = radial_offset(math.pi / 2.0, d)
        assert dx == pytest.approx(0.0, abs=1e-9)
        assert dy == pytest.approx(d, abs=1e-9)


# ---------------------------------------------------------------------------
# Tests: draw functions (smoke tests - no crash)
# ---------------------------------------------------------------------------

class TestDrawBranches:
    def test_draw_branches_no_error(self):
        tree = _make_tree(NEWICK_6)
        parent_map = _build_parent_map(tree)
        node_x = _build_node_x_phylogram(tree)
        coords = compute_circular_coords(tree, node_x, parent_map)

        fig, ax = plt.subplots()
        draw_circular_branches(ax, tree, coords, parent_map)
        plt.close(fig)

    def test_draw_tip_labels_no_error(self):
        tree = _make_tree(NEWICK_6)
        parent_map = _build_parent_map(tree)
        node_x = _build_node_x_phylogram(tree)
        coords = compute_circular_coords(tree, node_x, parent_map)

        fig, ax = plt.subplots()
        draw_circular_tip_labels(ax, tree, coords)
        plt.close(fig)

    def test_draw_colored_branch_no_error(self):
        pc = {"angle": 0.0, "radius": 1.0}
        cc = {"angle": 0.3, "radius": 2.0}
        fig, ax = plt.subplots()
        draw_circular_colored_branch(ax, pc, cc, "red")
        plt.close(fig)

    def test_draw_colored_arc_no_error(self):
        fig, ax = plt.subplots()
        draw_circular_colored_arc(ax, 0.0, 0.0, 1.0, 0.0, math.pi / 2, "blue")
        plt.close(fig)

    def test_draw_gradient_branch_no_error(self):
        import matplotlib.cm as cm
        pc = {"angle": 0.5, "radius": 1.0}
        cc = {"angle": 0.5, "radius": 3.0}
        fig, ax = plt.subplots()
        draw_circular_gradient_branch(ax, pc, cc, cm.viridis, 0.0, 1.0, 0.2, 0.8)
        plt.close(fig)

    def test_draw_multi_segment_branch_no_error(self):
        pc = {"angle": 1.0, "radius": 0.5}
        cc = {"angle": 1.0, "radius": 2.5}
        segments = [(0.0, 0.5, "A"), (0.5, 1.0, "B")]
        state_colors = {"A": "red", "B": "blue"}
        fig, ax = plt.subplots()
        draw_circular_multi_segment_branch(ax, pc, cc, segments, state_colors)
        plt.close(fig)
