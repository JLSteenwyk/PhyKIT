"""Tests for phykit.helpers.color_annotations."""

import os
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pytest
from Bio import Phylo
from io import StringIO

from phykit.helpers.color_annotations import (
    parse_color_file,
    resolve_mrca,
    get_clade_tip_ids,
    get_clade_branch_ids,
    draw_range_rect,
    draw_range_wedge,
)


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

    def test_parse_file_not_found(self):
        with pytest.raises(SystemExit):
            parse_color_file("/nonexistent/path/colors.tsv")


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


# ---------------------------------------------------------------------------
# draw_range_rect (smoke test)
# ---------------------------------------------------------------------------


class TestDrawRangeRect:
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
