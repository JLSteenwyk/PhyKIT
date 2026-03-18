import json
import os
import sys
import tempfile

import pytest
import numpy as np
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch

from phykit.services.tree.tree_space import TreeSpace
from phykit.errors import PhykitUserError

here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
GENE_TREES = str(SAMPLE_FILES / "gene_trees_simple.nwk")


def _base_args(**overrides):
    defaults = dict(
        trees=GENE_TREES,
        output=os.path.join(tempfile.mkdtemp(), "test_output.png"),
        metric="rf",
        method="mds",
        species_tree=None,
        k=None,
        seed=42,
        json=False,
        fig_width=None,
        fig_height=None,
        dpi=300,
        no_title=False,
        title=None,
        legend_position=None,
        ylabel_fontsize=None,
        xlabel_fontsize=None,
        title_fontsize=None,
        axis_fontsize=None,
        colors=None,
        ladderize=False,
        cladogram=False,
        circular=False,
        color_file=None,
    )
    defaults.update(overrides)
    return Namespace(**defaults)


@pytest.fixture
def default_args(tmp_path):
    return _base_args(output=str(tmp_path / "test_output.png"))


@pytest.fixture
def kf_args(tmp_path):
    return _base_args(output=str(tmp_path / "test_kf.png"), metric="kf")


@pytest.fixture
def species_tree_args(tmp_path):
    return _base_args(
        output=str(tmp_path / "test_sp.png"),
        species_tree=TREE_SIMPLE,
    )


class TestTreeParsing:
    def test_loads_all_gene_trees(self, default_args):
        svc = TreeSpace(default_args)
        trees = svc._parse_trees_from_source(GENE_TREES)
        assert len(trees) == 10

    def test_file_not_found(self, default_args):
        svc = TreeSpace(default_args)
        with pytest.raises(PhykitUserError):
            svc._parse_trees_from_source("no_such_file.nwk")

    def test_shared_taxa(self, default_args):
        svc = TreeSpace(default_args)
        trees = svc._parse_trees_from_source(GENE_TREES)
        shared = svc._get_shared_taxa(trees)
        assert len(shared) == 8


class TestDistanceMatrix:
    def test_rf_distance_matrix_symmetric(self, default_args):
        svc = TreeSpace(default_args)
        trees = svc._parse_trees_from_source(GENE_TREES)
        shared = svc._get_shared_taxa(trees)
        dist = svc._build_distance_matrix(trees, shared, "rf")
        np.testing.assert_array_almost_equal(dist, dist.T)

    def test_rf_distance_matrix_diagonal_zero(self, default_args):
        svc = TreeSpace(default_args)
        trees = svc._parse_trees_from_source(GENE_TREES)
        shared = svc._get_shared_taxa(trees)
        dist = svc._build_distance_matrix(trees, shared, "rf")
        np.testing.assert_array_almost_equal(np.diag(dist), 0.0)

    def test_rf_identical_trees_zero_distance(self, tmp_path, default_args):
        # Write two identical trees
        tree_file = tmp_path / "identical.nwk"
        tree_str = "((A,B),(C,D),(E,F));\n" * 5
        tree_file.write_text(tree_str)

        args = _base_args(
            trees=str(tree_file),
            output=str(tmp_path / "out.png"),
        )
        svc = TreeSpace(args)
        trees = svc._parse_trees_from_source(str(tree_file))
        shared = svc._get_shared_taxa(trees)
        dist = svc._build_distance_matrix(trees, shared, "rf")
        # All distances should be zero for identical trees
        np.testing.assert_array_almost_equal(dist, 0.0)

    def test_kf_distance_matrix(self, kf_args):
        svc = TreeSpace(kf_args)
        trees = svc._parse_trees_from_source(GENE_TREES)
        shared = svc._get_shared_taxa(trees)
        dist = svc._build_distance_matrix(trees, shared, "kf")
        # KF distance matrix should be symmetric
        np.testing.assert_array_almost_equal(dist, dist.T)
        # Diagonal should be zero
        np.testing.assert_array_almost_equal(np.diag(dist), 0.0)
        # KF distances should be non-negative
        assert np.all(dist >= 0)


class TestDimensionalityReduction:
    def test_mds_output_shape(self, default_args):
        svc = TreeSpace(default_args)
        trees = svc._parse_trees_from_source(GENE_TREES)
        shared = svc._get_shared_taxa(trees)
        dist = svc._build_distance_matrix(trees, shared, "rf")
        coords = svc._reduce_dimensions(dist, "mds", 42)
        assert coords.shape == (10, 2)

    def test_tsne_output_shape(self, default_args):
        svc = TreeSpace(default_args)
        trees = svc._parse_trees_from_source(GENE_TREES)
        shared = svc._get_shared_taxa(trees)
        dist = svc._build_distance_matrix(trees, shared, "rf")
        coords = svc._reduce_dimensions(dist, "tsne", 42)
        assert coords.shape == (10, 2)


class TestClustering:
    def test_clustering_labels_valid(self, default_args):
        svc = TreeSpace(default_args)
        trees = svc._parse_trees_from_source(GENE_TREES)
        shared = svc._get_shared_taxa(trees)
        dist = svc._build_distance_matrix(trees, shared, "rf")
        labels, k = svc._spectral_cluster(dist, None, 42)
        assert len(labels) == 10
        assert all(0 <= l < k for l in labels)
        assert k >= 2

    def test_auto_k_detection(self, default_args):
        svc = TreeSpace(default_args)
        trees = svc._parse_trees_from_source(GENE_TREES)
        shared = svc._get_shared_taxa(trees)
        dist = svc._build_distance_matrix(trees, shared, "rf")
        labels, k = svc._spectral_cluster(dist, None, 42)
        # k should be a reasonable number
        assert 2 <= k <= len(trees) // 2


class TestSpeciesTree:
    def test_species_tree_included(self, species_tree_args):
        svc = TreeSpace(species_tree_args)
        trees = svc._parse_trees_from_source(GENE_TREES)
        from Bio import Phylo
        sp_tree = Phylo.read(TREE_SIMPLE, "newick")
        shared = svc._get_shared_taxa(trees, sp_tree)
        all_trees = list(trees) + [sp_tree]
        dist = svc._build_distance_matrix(all_trees, shared, "rf")
        coords = svc._reduce_dimensions(dist, "mds", 42)
        # Should have 11 points (10 gene trees + 1 species tree)
        assert coords.shape == (11, 2)


class TestEndToEnd:
    @patch("builtins.print")
    def test_creates_png(self, mocked_print, tmp_path):
        output = tmp_path / "tree_space.png"
        args = _base_args(output=str(output), seed=42)
        svc = TreeSpace(args)
        svc.run()
        assert output.exists()
        assert output.stat().st_size > 0

    @patch("builtins.print")
    def test_json_output(self, mocked_print, tmp_path):
        output = tmp_path / "tree_space.png"
        args = _base_args(output=str(output), seed=42, json=True)
        svc = TreeSpace(args)
        svc.run()

        # Find the JSON output call
        json_str = mocked_print.call_args.args[0]
        payload = json.loads(json_str)
        assert payload["method"] == "mds"
        assert payload["metric"] == "rf"
        assert payload["n_trees"] == 10
        assert payload["n_taxa"] == 8
        assert "n_clusters" in payload
        assert "clusters" in payload
        assert "coordinates" in payload
        assert len(payload["coordinates"]) == 10
        assert payload["output_file"] == str(output)
