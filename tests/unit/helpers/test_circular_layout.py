"""Tests for phykit.helpers.circular_layout."""

import math
import subprocess
import sys
from io import StringIO

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest
from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin

import phykit.helpers.circular_layout as circular_layout_module
from phykit.helpers.circular_layout import (
    _preorder_clades_direct,
    _terminal_clades_direct,
    circular_branch_points,
    compute_circular_coords,
    draw_circular_branches,
    draw_circular_colored_arc,
    draw_circular_colored_arcs,
    draw_circular_colored_branch,
    draw_circular_gradient_branch,
    draw_circular_gradient_branches,
    draw_circular_multi_segment_branch,
    draw_circular_tip_labels,
    radial_offset,
)


def test_module_import_does_not_import_numpy():
    code = """
import sys
import phykit.helpers.circular_layout as module
assert hasattr(module.np, "__getattr__")
assert "typing" not in sys.modules
assert "numpy" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_lazy_numpy_caches_module_and_attributes():
    lazy_np = circular_layout_module._LazyNumpy()

    first_cos = lazy_np.cos
    second_cos = lazy_np.cos

    assert lazy_np._module is not None
    assert first_cos is second_cos
    assert lazy_np.__dict__["cos"] is first_cos


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

    def test_compute_circular_coords_uses_direct_tree_traversal(self, monkeypatch):
        tree = _make_tree("((A:1,B:1):1,C:2);")
        parent_map = _build_parent_map(tree)
        node_x = _build_node_x_phylogram(tree)
        tips = tree.get_terminals()
        internal = parent_map[id(tips[0])]

        def fail_traversal(*args, **kwargs):
            raise AssertionError("generic tree traversal should not be used")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)
        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)

        coords = compute_circular_coords(tree, node_x, parent_map)

        assert coords[id(tips[0])]["angle"] == pytest.approx(0.0, abs=1e-12)
        assert coords[id(tips[1])]["angle"] == pytest.approx(2.0 * math.pi / 3.0)
        assert coords[id(tips[2])]["angle"] == pytest.approx(4.0 * math.pi / 3.0)
        assert coords[id(internal)]["angle"] == pytest.approx(math.pi / 3.0)

    def test_precomputed_clade_lists_preserve_coordinates(self):
        tree = _make_tree(NEWICK_6)
        parent_map = _build_parent_map(tree)
        node_x = _build_node_x_phylogram(tree)
        preorder_clades = list(tree.find_clades(order="preorder"))
        terminal_clades = tree.get_terminals()

        expected = compute_circular_coords(tree, node_x, parent_map)
        observed = compute_circular_coords(
            tree,
            node_x,
            parent_map,
            preorder_clades=preorder_clades,
            terminal_clades=terminal_clades,
        )

        assert set(observed) == set(expected)
        for clade_id, expected_coords in expected.items():
            for key, expected_value in expected_coords.items():
                assert observed[clade_id][key] == pytest.approx(expected_value)


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

    def test_branch_points_compute_trig_once(self, monkeypatch):
        cos_calls = 0
        sin_calls = 0
        real_cos = math.cos
        real_sin = math.sin

        def track_cos(angle):
            nonlocal cos_calls
            cos_calls += 1
            return real_cos(angle)

        def track_sin(angle):
            nonlocal sin_calls
            sin_calls += 1
            return real_sin(angle)

        monkeypatch.setattr(math, "cos", track_cos)
        monkeypatch.setattr(math, "sin", track_sin)

        pts = circular_branch_points(
            {"angle": 0.5, "radius": 1.0},
            {"angle": 0.5, "radius": 2.0},
            20,
        )

        assert len(pts) == 20
        assert cos_calls == 1
        assert sin_calls == 1


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

    def test_draw_branches_batches_real_matplotlib_axes(self, monkeypatch):
        tree = _make_tree(NEWICK_6)
        parent_map = _build_parent_map(tree)
        node_x = _build_node_x_phylogram(tree)
        coords = compute_circular_coords(tree, node_x, parent_map)

        fig, ax = plt.subplots()

        def fail_plot(*_args, **_kwargs):
            raise AssertionError("real axes should use line collections")

        monkeypatch.setattr(ax, "plot", fail_plot)

        draw_circular_branches(ax, tree, coords, parent_map)

        assert len(ax.collections) == 2
        plt.close(fig)

    def test_draw_branches_uses_direct_traversal(self, monkeypatch):
        tree = _make_tree(NEWICK_6)
        parent_map = _build_parent_map(tree)
        node_x = _build_node_x_phylogram(tree)
        coords = compute_circular_coords(tree, node_x, parent_map)

        def fail_find_clades(*args, **kwargs):
            raise AssertionError("generic tree traversal should not be used")

        class NoopAxes:
            def __init__(self):
                self.plot_calls = 0

            def plot(self, *args, **kwargs):
                self.plot_calls += 1

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        ax = NoopAxes()
        draw_circular_branches(ax, tree, coords, parent_map)

        internal_branch_count = 0
        stack = [tree.root]
        while stack:
            clade = stack.pop()
            children = clade.clades
            if len(children) >= 2:
                internal_branch_count += 1
            child_count = len(children)
            if child_count == 2:
                stack.append(children[1])
                stack.append(children[0])
            else:
                for index in range(child_count - 1, -1, -1):
                    stack.append(children[index])

        assert ax.plot_calls == len(parent_map) + internal_branch_count

    def test_draw_tip_labels_no_error(self):
        tree = _make_tree(NEWICK_6)
        parent_map = _build_parent_map(tree)
        node_x = _build_node_x_phylogram(tree)
        coords = compute_circular_coords(tree, node_x, parent_map)

        fig, ax = plt.subplots()
        draw_circular_tip_labels(ax, tree, coords)
        plt.close(fig)

    def test_draw_tip_labels_emits_expected_text_calls(self):
        tree = _make_tree("(A:1,B:1,C:1,D:1);")
        tips = tree.get_terminals()
        coords = {
            id(tips[0]): {"angle": 0.0, "radius": 1.0},
            id(tips[1]): {"angle": math.pi, "radius": 1.0},
            id(tips[2]): {"angle": math.pi / 2, "radius": 1.0},
            id(tips[3]): {"angle": -math.pi / 2, "radius": 1.0},
        }

        class TextAxes:
            def __init__(self):
                self.calls = []

            def text(self, *args, **kwargs):
                self.calls.append((args, kwargs))

        ax = TextAxes()
        draw_circular_tip_labels(ax, tree, coords, fontsize=7, offset=0.5)

        assert [call[0][2] for call in ax.calls] == ["A", "B", "C", "D"]
        assert [call[1]["ha"] for call in ax.calls] == [
            "left",
            "right",
            "left",
            "left",
        ]
        assert [call[1]["rotation"] for call in ax.calls] == [
            0.0,
            0.0,
            90.0,
            -90.0,
        ]
        assert all(call[1]["fontsize"] == 7 for call in ax.calls)
        assert all(call[1]["rotation_mode"] == "anchor" for call in ax.calls)

    def test_draw_tip_labels_uses_direct_terminal_traversal(self, monkeypatch):
        tree = _make_tree(NEWICK_6)
        parent_map = _build_parent_map(tree)
        node_x = _build_node_x_phylogram(tree)
        coords = compute_circular_coords(tree, node_x, parent_map)

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("standard trees should not call get_terminals")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        fig, ax = plt.subplots()
        draw_circular_tip_labels(ax, tree, coords)
        plt.close(fig)

    def test_direct_clade_helpers_preserve_mixed_child_order(self, monkeypatch):
        tree = _make_tree("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1,G:1);")
        expected_preorder = list(tree.find_clades(order="preorder"))
        expected_terminals = tree.get_terminals()

        def fail_traversal(*args, **kwargs):
            raise AssertionError("direct clade helpers should not use generic traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)

        assert _preorder_clades_direct(tree) == expected_preorder
        assert _terminal_clades_direct(tree) == expected_terminals

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

    def test_draw_colored_arc_uses_vectorized_points(self, monkeypatch):
        class NoopAxes:
            def __init__(self):
                self.args = None

            def plot(self, *args, **kwargs):
                self.args = args

        def fail_scalar_trig(*_args, **_kwargs):
            raise AssertionError("arc point generation should use NumPy trig")

        monkeypatch.setattr(math, "cos", fail_scalar_trig)
        monkeypatch.setattr(math, "sin", fail_scalar_trig)

        ax = NoopAxes()
        draw_circular_colored_arc(ax, 0.0, 0.0, 1.0, 0.0, math.pi / 2, "blue")

        xs, ys = ax.args
        assert isinstance(xs, np.ndarray)
        assert isinstance(ys, np.ndarray)
        assert xs.shape == (61,)
        assert ys.shape == (61,)

    def test_draw_colored_arcs_batches_real_axes(self, monkeypatch):
        import matplotlib.axes
        from matplotlib.collections import LineCollection

        arcs = [
            (0.0, 0.0, 1.0, 0.0, math.pi / 2, "blue"),
            (0.0, 0.0, 2.0, math.pi / 2, math.pi, "red"),
        ]
        fig, ax = plt.subplots()
        original_add_collection = matplotlib.axes.Axes.add_collection
        line_collections = []

        def fail_plot(*args, **kwargs):
            raise AssertionError("colored arcs should be batched")

        def count_collection(self, collection, *args, **kwargs):
            if isinstance(collection, LineCollection):
                line_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "plot", fail_plot)
        monkeypatch.setattr(
            matplotlib.axes.Axes, "add_collection", count_collection
        )

        draw_circular_colored_arcs(ax, arcs)

        assert len(line_collections) == 1
        assert len(line_collections[0].get_segments()) == 2
        plt.close(fig)

    def test_draw_gradient_branch_no_error(self):
        import matplotlib.cm as cm
        pc = {"angle": 0.5, "radius": 1.0}
        cc = {"angle": 0.5, "radius": 3.0}
        fig, ax = plt.subplots()
        draw_circular_gradient_branch(ax, pc, cc, cm.viridis, 0.0, 1.0, 0.2, 0.8)
        plt.close(fig)

    def test_draw_gradient_branch_batches_real_axes(self, monkeypatch):
        import matplotlib.axes
        import matplotlib.cm as cm
        from matplotlib.collections import LineCollection

        pc = {"angle": 0.5, "radius": 1.0}
        cc = {"angle": 0.5, "radius": 3.0}
        fig, ax = plt.subplots()
        original_add_collection = matplotlib.axes.Axes.add_collection
        line_collections = []

        def fail_plot(*args, **kwargs):
            raise AssertionError("gradient branch segments should be batched")

        def count_collection(self, collection, *args, **kwargs):
            if isinstance(collection, LineCollection):
                line_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "plot", fail_plot)
        monkeypatch.setattr(
            matplotlib.axes.Axes, "add_collection", count_collection
        )

        draw_circular_gradient_branch(
            ax, pc, cc, cm.viridis, 0.0, 1.0, 0.2, 0.8
        )

        assert len(line_collections) == 1
        plt.close(fig)

    def test_draw_gradient_branches_batches_all_branches_real_axes(self, monkeypatch):
        import matplotlib.axes
        import matplotlib.cm as cm
        from matplotlib.collections import LineCollection

        branch_data = [
            (
                {"angle": 0.5, "radius": 1.0},
                {"angle": 0.5, "radius": 3.0},
                0.2,
                0.8,
            ),
            (
                {"angle": 1.0, "radius": 0.5},
                {"angle": 1.0, "radius": 2.0},
                0.8,
                0.1,
            ),
        ]
        fig, ax = plt.subplots()
        original_add_collection = matplotlib.axes.Axes.add_collection
        line_collections = []

        def fail_plot(*args, **kwargs):
            raise AssertionError("gradient branches should be batched together")

        def count_collection(self, collection, *args, **kwargs):
            if isinstance(collection, LineCollection):
                line_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "plot", fail_plot)
        monkeypatch.setattr(
            matplotlib.axes.Axes, "add_collection", count_collection
        )

        draw_circular_gradient_branches(
            ax, branch_data, cm.viridis, 0.0, 1.0
        )

        assert len(line_collections) == 1
        assert len(line_collections[0].get_segments()) == 60
        plt.close(fig)

    def test_draw_multi_segment_branch_no_error(self):
        pc = {"angle": 1.0, "radius": 0.5}
        cc = {"angle": 1.0, "radius": 2.5}
        segments = [(0.0, 0.5, "A"), (0.5, 1.0, "B")]
        state_colors = {"A": "red", "B": "blue"}
        fig, ax = plt.subplots()
        draw_circular_multi_segment_branch(ax, pc, cc, segments, state_colors)
        plt.close(fig)

    def test_draw_multi_segment_branch_batches_real_axes(self, monkeypatch):
        import matplotlib.axes
        from matplotlib.collections import LineCollection

        pc = {"angle": 1.0, "radius": 0.5}
        cc = {"angle": 1.0, "radius": 2.5}
        segments = [(0.0, 0.25, "A"), (0.25, 0.75, "B"), (0.75, 1.0, "A")]
        state_colors = {"A": "red", "B": "blue"}
        fig, ax = plt.subplots()
        original_add_collection = matplotlib.axes.Axes.add_collection
        line_collections = []

        def fail_plot(*args, **kwargs):
            raise AssertionError("multi-segment branches should be batched")

        def count_collection(self, collection, *args, **kwargs):
            if isinstance(collection, LineCollection):
                line_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "plot", fail_plot)
        monkeypatch.setattr(
            matplotlib.axes.Axes, "add_collection", count_collection
        )

        draw_circular_multi_segment_branch(ax, pc, cc, segments, state_colors)

        assert len(line_collections) == 1
        plt.close(fig)
