"""Tests for phykit.helpers.color_annotations."""

import os
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pytest
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, TreeMixin
from io import StringIO

from phykit.helpers.color_annotations import (
    _range_wedge_angle_bounds,
    _range_rect_tip_bounds,
    _terminal_clades,
    _valid_mrca_taxa,
    apply_label_colors,
    parse_color_file,
    resolve_mrca,
    get_clade_tip_ids,
    get_clade_branch_ids,
    draw_range_rect,
    draw_range_wedge,
)


def _legacy_range_wedge_angle_bounds(sorted_angles):
    import math

    if len(sorted_angles) >= 2:
        gaps = [
            sorted_angles[i + 1] - sorted_angles[i]
            for i in range(len(sorted_angles) - 1)
        ]
        min_gap = min(gaps) if gaps else 0.05
        pad = min_gap * 0.5
    else:
        pad = 0.05

    complement_gap = (2 * math.pi) - (sorted_angles[-1] - sorted_angles[0])
    all_gaps = [
        sorted_angles[i + 1] - sorted_angles[i]
        for i in range(len(sorted_angles) - 1)
    ]
    all_gaps.append(complement_gap)

    max_gap_idx = all_gaps.index(max(all_gaps))
    if max_gap_idx < len(sorted_angles) - 1:
        angle_min = sorted_angles[max_gap_idx + 1] - pad
        angle_max = sorted_angles[max_gap_idx] + pad
        if angle_max < angle_min:
            angle_max += 2 * math.pi
    else:
        angle_min = sorted_angles[0] - pad
        angle_max = sorted_angles[-1] + pad
    return angle_min, angle_max


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_tree(newick):
    return Phylo.read(StringIO(newick), "newick")


def _write_tmp(content: str) -> str:
    fd, path = tempfile.mkstemp(suffix=".tsv")
    with os.fdopen(fd, "w") as fh:
        fh.write(content)
    return path


# ---------------------------------------------------------------------------
# parse_color_file
# ---------------------------------------------------------------------------


class TestParseColorFile:
    def test_parse_label_entry(self):
        path = _write_tmp("A\tlabel\t#ff0000\n")
        try:
            result = parse_color_file(path)
            assert result["labels"] == {"A": "#ff0000"}
            assert result["ranges"] == []
            assert result["clades"] == []
        finally:
            os.unlink(path)

    def test_parse_range_entry(self):
        path = _write_tmp("A,B\trange\t#00ff00\tMyGroup\n")
        try:
            result = parse_color_file(path)
            assert result["labels"] == {}
            assert len(result["ranges"]) == 1
            taxa, color, label = result["ranges"][0]
            assert taxa == ["A", "B"]
            assert color == "#00ff00"
            assert label == "MyGroup"
        finally:
            os.unlink(path)

    def test_parse_clade_entry(self):
        path = _write_tmp("A,B\tclade\t#0000ff\n")
        try:
            result = parse_color_file(path)
            assert result["labels"] == {}
            assert result["ranges"] == []
            assert len(result["clades"]) == 1
            taxa, color, label = result["clades"][0]
            assert taxa == ["A", "B"]
            assert color == "#0000ff"
            assert label is None
        finally:
            os.unlink(path)

    def test_parse_comments_skipped(self):
        path = _write_tmp("# This is a comment\nA\tlabel\tred\n")
        try:
            result = parse_color_file(path)
            assert result["labels"] == {"A": "red"}
        finally:
            os.unlink(path)

    def test_parse_empty_lines_skipped(self):
        path = _write_tmp("\n\nA\tlabel\tblue\n\n")
        try:
            result = parse_color_file(path)
            assert result["labels"] == {"A": "blue"}
        finally:
            os.unlink(path)

    def test_parse_missing_label_field(self):
        path = _write_tmp("A,B\trange\t#00ff00\n")
        try:
            result = parse_color_file(path)
            assert len(result["ranges"]) == 1
            _, _, label = result["ranges"][0]
            assert label is None
        finally:
            os.unlink(path)

    def test_parse_mixed_case_types_and_ignores_extra_columns(self):
        path = _write_tmp(
            "A\tLabel\tred\tignored\n"
            "B,C\tRange\tblue\tRange label\textra\n"
            "D,E\tClade\tgreen\tClade label\textra\n"
        )
        try:
            result = parse_color_file(path)
            assert result["labels"] == {"A": "red"}
            assert result["ranges"] == [(["B", "C"], "blue", "Range label")]
            assert result["clades"] == [(["D", "E"], "green", "Clade label")]
        finally:
            os.unlink(path)

    def test_parse_consolidated_type_dispatch_ignores_unknown_types(self):
        path = _write_tmp(
            "A\tlabel\tred\n"
            "B,C\tRange\tblue\tRange label\n"
            "D,E\tclade\tgreen\tClade label\n"
            "F\tunknown\tblack\tignored\n"
        )
        try:
            result = parse_color_file(path)
            assert result["labels"] == {"A": "red"}
            assert result["ranges"] == [(["B", "C"], "blue", "Range label")]
            assert result["clades"] == [(["D", "E"], "green", "Clade label")]
        finally:
            os.unlink(path)

    def test_parse_color_file_streams_without_readlines(self, monkeypatch):
        class StreamingOnlyFile:
            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, traceback):
                return False

            def __iter__(self):
                return iter(
                    [
                        "# ignored\n",
                        "\n",
                        "A\tlabel\tred\n",
                        "A,B\trange\tblue\tRange label\n",
                        "C,D\tclade\tgreen\n",
                    ]
                )

            def readlines(self):
                raise AssertionError("parse_color_file should stream rows")

        monkeypatch.setattr("builtins.open", lambda path: StreamingOnlyFile())

        result = parse_color_file("colors.tsv")

        assert result["labels"] == {"A": "red"}
        assert result["ranges"] == [(["A", "B"], "blue", "Range label")]
        assert result["clades"] == [(["C", "D"], "green", None)]

    def test_parse_file_not_found(self):
        with pytest.raises(SystemExit):
            parse_color_file("/nonexistent/path/colors.tsv")


class TestApplyLabelColors:
    def test_applies_matching_label_colors(self):
        fig, ax = plt.subplots()
        try:
            ax.text(0, 0, "A")
            ax.text(0, 1, "B")
            ax.text(0, 2, "C")

            count = apply_label_colors(ax, {"A": "#ff0000", "C": "#0000ff"})

            assert count == 2
            assert ax.texts[0].get_color() == "#ff0000"
            assert ax.texts[1].get_color() == "black"
            assert ax.texts[2].get_color() == "#0000ff"
        finally:
            plt.close(fig)

    def test_scans_text_artists_once(self, monkeypatch):
        fig, ax = plt.subplots()
        try:
            for index in range(6):
                ax.text(0, index, f"T{index}")

            get_text_calls = 0
            original_get_text = ax.texts[0].__class__.get_text

            def count_get_text(self):
                nonlocal get_text_calls
                get_text_calls += 1
                return original_get_text(self)

            monkeypatch.setattr(ax.texts[0].__class__, "get_text", count_get_text)

            count = apply_label_colors(
                ax,
                {f"T{index}": f"C{index}" for index in range(6)},
            )

            assert count == 6
            assert get_text_calls == 6
        finally:
            plt.close(fig)


# ---------------------------------------------------------------------------
# resolve_mrca
# ---------------------------------------------------------------------------


class TestResolveMrca:
    def test_resolve_mrca_two_taxa(self):
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        mrca = resolve_mrca(tree, ["A", "B"])
        assert mrca is not None
        tip_names = sorted(t.name for t in mrca.get_terminals())
        assert tip_names == ["A", "B"]

    def test_resolve_mrca_missing_taxon(self):
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        # "X" is not in tree, but A and B are still valid
        mrca = resolve_mrca(tree, ["A", "B", "X"])
        assert mrca is not None
        tip_names = sorted(t.name for t in mrca.get_terminals())
        assert tip_names == ["A", "B"]

    def test_resolve_mrca_single_valid_taxon(self):
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        mrca = resolve_mrca(tree, ["A", "X", "Y"])
        assert mrca is None

    def test_resolve_mrca_uses_direct_tip_name_traversal(self, monkeypatch):
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")

        def fail_get_terminals(self):
            raise AssertionError("get_terminals should not be used")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        mrca = resolve_mrca(tree, ["A", "B"])
        assert mrca is not None

    def test_resolve_mrca_stops_after_all_requested_taxa_are_found(self):
        class ExplodingClade:
            @property
            def clades(self):
                raise AssertionError("unneeded clade should not be traversed")

        class SimpleClade:
            def __init__(self, name=None, clades=None):
                self.name = name
                self.clades = clades or []

        class SimpleTree:
            def __init__(self):
                self.root = SimpleClade(
                    clades=[
                        SimpleClade(clades=[SimpleClade("A"), SimpleClade("B")]),
                        ExplodingClade(),
                    ]
                )
                self.requested_taxa = None
                self.mrca = object()

            def common_ancestor(self, taxa):
                self.requested_taxa = taxa
                return self.mrca

        tree = SimpleTree()

        assert resolve_mrca(tree, ["A", "B"]) is tree.mrca
        assert tree.requested_taxa == ["A", "B"]

    def test_valid_mrca_taxa_handles_mixed_child_counts(self, monkeypatch):
        tree = _make_tree("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1,G:1);")

        def fail_get_terminals(self):
            raise AssertionError("standard clades should not use get_terminals")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        valid = _valid_mrca_taxa(tree.root, ["F", "A", "missing", "C"])

        assert valid == ["F", "A", "C"]


# ---------------------------------------------------------------------------
# get_clade_tip_ids
# ---------------------------------------------------------------------------


class TestGetCladeTipIds:
    def test_get_clade_tip_ids(self):
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        mrca = tree.common_ancestor(["A", "B"])
        tip_ids = get_clade_tip_ids(mrca)
        # Should contain exactly 2 tip ids
        assert len(tip_ids) == 2
        actual_ids = {id(t) for t in mrca.get_terminals()}
        assert tip_ids == actual_ids

    def test_get_clade_tip_ids_uses_direct_traversal(self, monkeypatch):
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        mrca = tree.common_ancestor(["A", "B"])
        expected = {
            id(child)
            for child in mrca.clades
            if not child.clades
        }

        def fail_get_terminals(self):
            raise AssertionError("get_terminals should not be used")

        monkeypatch.setattr(Clade, "get_terminals", fail_get_terminals)

        assert get_clade_tip_ids(mrca) == expected

    def test_terminal_clades_preserves_mixed_child_order(self, monkeypatch):
        tree = _make_tree("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1,G:1);")

        def fail_get_terminals(self):
            raise AssertionError("standard clades should not use get_terminals")

        monkeypatch.setattr(Clade, "get_terminals", fail_get_terminals)

        assert [tip.name for tip in _terminal_clades(tree.root)] == [
            "A",
            "B",
            "C",
            "D",
            "E",
            "F",
            "G",
        ]


# ---------------------------------------------------------------------------
# get_clade_branch_ids
# ---------------------------------------------------------------------------


class TestGetCladeBranchIds:
    def test_get_clade_branch_ids(self):
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        mrca = tree.common_ancestor(["A", "B"])
        # Build parent_map
        parent_map = {}
        for clade in tree.find_clades(order="preorder"):
            for child in clade.clades:
                parent_map[id(child)] = clade

        branch_ids = get_clade_branch_ids(tree, mrca, parent_map)
        # Should include the mrca itself plus A and B = 3 nodes
        assert id(mrca) in branch_ids
        assert len(branch_ids) == 3

    def test_get_clade_branch_ids_uses_direct_traversal(self, monkeypatch):
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        mrca = tree.common_ancestor(["A", "B"])

        def fail_find_clades(self, *args, **kwargs):
            raise AssertionError("standard clade branch IDs should not use find_clades")

        monkeypatch.setattr(Clade, "find_clades", fail_find_clades)

        branch_ids = get_clade_branch_ids(tree, mrca, {})

        assert branch_ids == {
            id(mrca),
            id(mrca.clades[0]),
            id(mrca.clades[1]),
        }

    def test_get_clade_branch_ids_handles_mixed_child_counts(self, monkeypatch):
        tree = _make_tree("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);")

        def fail_find_clades(self, *args, **kwargs):
            raise AssertionError("standard clade branch IDs should not use find_clades")

        monkeypatch.setattr(Clade, "find_clades", fail_find_clades)

        branch_ids = get_clade_branch_ids(tree, tree.root, {})

        expected = {id(tree.root)}
        for child in tree.root.clades:
            expected.add(id(child))
            for grandchild in child.clades:
                expected.add(id(grandchild))

        assert branch_ids == expected


# ---------------------------------------------------------------------------
# draw_range_rect (smoke test)
# ---------------------------------------------------------------------------


class TestDrawRangeRect:
    def test_range_rect_tip_bounds_scans_tips_once_and_keeps_axes_independent(self):
        class Tip:
            pass

        class SinglePassTips:
            def __init__(self, tips):
                self.tips = tips
                self.iterations = 0

            def __iter__(self):
                self.iterations += 1
                if self.iterations > 1:
                    raise AssertionError("tips should be scanned once")
                return iter(self.tips)

        tip_a = Tip()
        tip_b = Tip()
        tip_c = Tip()
        tips = SinglePassTips([tip_a, tip_b, tip_c])
        node_x = {id(tip_a): 3.0, id(tip_c): 1.0}
        node_y = {id(tip_b): 4.0, id(tip_c): 2.0}

        assert _range_rect_tip_bounds(tips, node_x, node_y) == (
            1.0,
            3.0,
            2.0,
            4.0,
        )
        assert tips.iterations == 1
        assert _range_rect_tip_bounds([], node_x, node_y) is None
        assert _range_rect_tip_bounds([tip_a], node_x, {}) is None

    def test_draw_range_rect_no_error(self):
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        mrca = tree.common_ancestor(["A", "B"])

        # Build node_x and node_y dicts
        node_x = {}
        node_y = {}
        y_counter = 0
        for clade in tree.find_clades(order="preorder"):
            if clade.is_terminal():
                node_x[id(clade)] = 2.0
                node_y[id(clade)] = float(y_counter)
                y_counter += 1
            else:
                node_x[id(clade)] = 0.5
                node_y[id(clade)] = float(y_counter) * 0.5

        fig, ax = plt.subplots()
        # Should not raise
        draw_range_rect(ax, tree, mrca, "#ff0000", node_x, node_y)
        plt.close(fig)


# ---------------------------------------------------------------------------
# draw_range_wedge (smoke test)
# ---------------------------------------------------------------------------


class TestDrawRangeWedge:
    def test_range_wedge_angle_bounds_matches_legacy_two_pass_calculation(self):
        cases = [
            [0.1],
            [0.1, 0.2, 0.35, 0.6],
            [0.05, 0.12, 6.05, 6.15],
            [0.2, 2.2],
            [idx * 0.01 for idx in range(100)],
        ]
        for angles in cases:
            assert _range_wedge_angle_bounds(angles) == pytest.approx(
                _legacy_range_wedge_angle_bounds(angles)
            )

    def test_draw_range_wedge_no_error(self):
        import math

        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        mrca = tree.common_ancestor(["A", "B"])

        # Build circular coords dict: id(node) -> {"x", "y", "angle", "radius"}
        coords = {}
        tips = list(tree.get_terminals())
        n = len(tips)
        for i, tip in enumerate(tips):
            angle = 2 * math.pi * i / n
            r = 2.0
            coords[id(tip)] = {
                "x": r * math.cos(angle),
                "y": r * math.sin(angle),
                "angle": angle,
                "radius": r,
            }

        # Add MRCA coord at smaller radius
        mrca_angle = math.pi / 4
        coords[id(mrca)] = {
            "x": 0.5 * math.cos(mrca_angle),
            "y": 0.5 * math.sin(mrca_angle),
            "angle": mrca_angle,
            "radius": 0.5,
        }

        fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
        # draw_range_wedge works on a regular (non-polar) axes too
        fig2, ax2 = plt.subplots()

        # Should not raise
        draw_range_wedge(ax2, tree, mrca, "#00ff00", coords)
        plt.close(fig)
        plt.close(fig2)

    def test_draw_range_wedge_uses_mrca_and_tip_radii(self):
        import math

        tree = _make_tree("((A:1,B:1):1,C:1);")
        mrca = tree.common_ancestor(["A", "B"])
        tip_a, tip_b = mrca.clades
        coords = {
            id(tip_a): {"angle": 0.0, "radius": 2.0},
            id(tip_b): {"angle": math.pi / 2, "radius": 3.0},
            id(mrca): {"angle": math.pi / 4, "radius": 0.5},
        }

        fig, ax = plt.subplots()
        draw_range_wedge(ax, tree, mrca, "#00ff00", coords)

        wedge = ax.patches[0]
        assert wedge.r == pytest.approx(3.125)
        assert wedge.width == pytest.approx(2.625)
        plt.close(fig)
