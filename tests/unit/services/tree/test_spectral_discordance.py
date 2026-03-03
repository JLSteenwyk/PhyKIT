import json
import os
import sys
import tempfile

import pytest
import numpy as np
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch

from phykit.services.tree.spectral_discordance import SpectralDiscordance
from phykit.errors import PhykitUserError

here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
GENE_TREES = str(SAMPLE_FILES / "gene_trees_simple.nwk")


@pytest.fixture
def default_args():
    return Namespace(
        tree=TREE_SIMPLE,
        gene_trees=GENE_TREES,
        metric="nrf",
        clusters=None,
        n_pcs=None,
        top_loadings=5,
        plot=None,
        json=False,
    )


@pytest.fixture
def no_species_tree_args():
    return Namespace(
        tree=None,
        gene_trees=GENE_TREES,
        metric="nrf",
        clusters=None,
        n_pcs=None,
        top_loadings=5,
        plot=None,
        json=False,
    )


@pytest.fixture
def wrf_args():
    return Namespace(
        tree=TREE_SIMPLE,
        gene_trees=GENE_TREES,
        metric="wrf",
        clusters=None,
        n_pcs=None,
        top_loadings=5,
        plot=None,
        json=False,
    )


class TestGeneTreeLoading:
    def test_loads_all_gene_trees(self, default_args):
        svc = SpectralDiscordance(default_args)
        trees = svc._parse_gene_trees(GENE_TREES)
        assert len(trees) == 10

    def test_file_not_found(self, default_args):
        svc = SpectralDiscordance(default_args)
        with pytest.raises(PhykitUserError):
            svc._parse_gene_trees("no_such_file.nwk")

    def test_shared_taxa(self, default_args):
        svc = SpectralDiscordance(default_args)
        trees = svc._parse_gene_trees(GENE_TREES)
        shared = svc._get_shared_taxa(trees, species_tree=None)
        assert len(shared) == 8


class TestBipartitionExtraction:
    def test_canonical_split_smaller_side(self, default_args):
        svc = SpectralDiscordance(default_args)
        all_taxa = frozenset({"A", "B", "C", "D", "E"})
        result = svc._canonical_split(frozenset({"A", "B"}), all_taxa)
        assert result == frozenset({"A", "B"})

    def test_canonical_split_larger_side_flips(self, default_args):
        svc = SpectralDiscordance(default_args)
        all_taxa = frozenset({"A", "B", "C", "D", "E"})
        result = svc._canonical_split(frozenset({"A", "B", "C", "D"}), all_taxa)
        assert result == frozenset({"E"})

    def test_extract_splits_returns_frozensets(self, default_args):
        svc = SpectralDiscordance(default_args)
        trees = svc._parse_gene_trees(GENE_TREES)
        shared = svc._get_shared_taxa(trees, species_tree=None)
        all_taxa_fs = frozenset(shared)
        splits = svc._extract_splits(trees[0], all_taxa_fs)
        assert isinstance(splits, set)
        for s in splits:
            assert isinstance(s, frozenset)

    def test_extract_splits_no_trivial(self, default_args):
        svc = SpectralDiscordance(default_args)
        trees = svc._parse_gene_trees(GENE_TREES)
        shared = svc._get_shared_taxa(trees, species_tree=None)
        all_taxa_fs = frozenset(shared)
        splits = svc._extract_splits(trees[0], all_taxa_fs)
        for s in splits:
            assert len(s) > 1 or len(all_taxa_fs) - len(s) > 1
            assert s != all_taxa_fs

    def test_bipartition_matrix_shape(self, default_args):
        svc = SpectralDiscordance(default_args)
        trees = svc._parse_gene_trees(GENE_TREES)
        shared = svc._get_shared_taxa(trees, species_tree=None)
        X, bip_index, _ = svc._build_bipartition_matrix(trees, shared, metric="nrf")
        assert X.shape[0] == 10
        assert X.shape[1] == len(bip_index)
        assert X.shape[1] > 0

    def test_bipartition_matrix_binary_for_nrf(self, default_args):
        svc = SpectralDiscordance(default_args)
        trees = svc._parse_gene_trees(GENE_TREES)
        shared = svc._get_shared_taxa(trees, species_tree=None)
        X, _, _ = svc._build_bipartition_matrix(trees, shared, metric="nrf")
        assert set(np.unique(X)).issubset({0.0, 1.0})

    def test_bipartition_matrix_wrf_has_branch_lengths(self, default_args):
        svc = SpectralDiscordance(default_args)
        trees = svc._parse_gene_trees(GENE_TREES)
        shared = svc._get_shared_taxa(trees, species_tree=None)
        X, _, _ = svc._build_bipartition_matrix(trees, shared, metric="wrf")
        assert np.max(X) > 1.0

    def test_species_tree_flags(self, default_args):
        svc = SpectralDiscordance(default_args)
        species_tree = svc.read_tree_file()
        trees = svc._parse_gene_trees(GENE_TREES)
        shared = svc._get_shared_taxa(trees, species_tree=species_tree)
        _, bip_index, sp_flags = svc._build_bipartition_matrix(
            trees, shared, metric="nrf", species_tree=species_tree
        )
        assert any(sp_flags.values())
        assert not all(sp_flags.values())


class TestPCA:
    def _get_pca(self, args):
        svc = SpectralDiscordance(args)
        trees = svc._parse_gene_trees(GENE_TREES)
        species_tree = svc.read_tree_file() if svc.tree_file_path else None
        shared = svc._get_shared_taxa(trees, species_tree)
        X, bip_index, sp_flags = svc._build_bipartition_matrix(
            trees, shared, metric=args.metric, species_tree=species_tree
        )
        scores, var_explained, loadings = svc._run_pca(X)
        return scores, var_explained, loadings, bip_index

    def test_scores_shape(self, default_args):
        scores, _, _, _ = self._get_pca(default_args)
        assert scores.shape[0] == 10
        assert scores.shape[1] <= 10

    def test_variance_explained_sums_to_one(self, default_args):
        _, ve, _, _ = self._get_pca(default_args)
        assert pytest.approx(np.sum(ve), abs=1e-6) == 1.0

    def test_variance_explained_descending(self, default_args):
        _, ve, _, _ = self._get_pca(default_args)
        for i in range(len(ve) - 1):
            assert ve[i] >= ve[i + 1] - 1e-10

    def test_loadings_shape(self, default_args):
        scores, _, loadings, bip_index = self._get_pca(default_args)
        assert loadings.shape[1] == len(bip_index)
        assert loadings.shape[0] == scores.shape[1]

    def test_first_pc_captures_most_variance(self, default_args):
        _, ve, _, _ = self._get_pca(default_args)
        assert ve[0] > 0.3

    def test_pca_wrf(self, wrf_args):
        scores, ve, loadings, _ = self._get_pca(wrf_args)
        assert scores.shape[0] == 10
        assert pytest.approx(np.sum(ve), abs=1e-6) == 1.0


class TestSpectralClustering:
    def _get_clusters(self, args, n_clusters=None):
        svc = SpectralDiscordance(args)
        trees = svc._parse_gene_trees(GENE_TREES)
        species_tree = svc.read_tree_file() if svc.tree_file_path else None
        shared = svc._get_shared_taxa(trees, species_tree)
        X, _, _ = svc._build_bipartition_matrix(trees, shared, metric=args.metric)
        X_centered = X - X.mean(axis=0)
        labels, K, eigengaps = svc._spectral_cluster(X_centered, n_clusters=n_clusters)
        return labels, K, eigengaps

    def test_labels_shape(self, default_args):
        labels, _, _ = self._get_clusters(default_args)
        assert len(labels) == 10

    def test_labels_valid_range(self, default_args):
        labels, K, _ = self._get_clusters(default_args)
        assert all(0 <= l < K for l in labels)

    def test_k_at_least_2(self, default_args):
        _, K, _ = self._get_clusters(default_args)
        assert K >= 2

    def test_eigengaps_positive(self, default_args):
        _, _, eigengaps = self._get_clusters(default_args)
        assert all(g >= -1e-10 for g in eigengaps)

    def test_override_k(self, default_args):
        labels, K, _ = self._get_clusters(default_args, n_clusters=3)
        assert K == 3
        assert len(set(labels)) <= 3

    def test_concordant_trees_cluster_together(self, default_args):
        labels, K, _ = self._get_clusters(default_args)
        assert K >= 2
        concordant_labels = [labels[i] for i in range(6)]
        from collections import Counter
        most_common = Counter(concordant_labels).most_common(1)[0][1]
        assert most_common >= 4


class TestRun:
    @patch("builtins.print")
    def test_text_output(self, mocked_print, default_args):
        svc = SpectralDiscordance(default_args)
        svc.run()
        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Spectral Discordance" in all_output
        assert "Variance explained" in all_output
        assert "cluster" in all_output.lower()

    @patch("builtins.print")
    def test_json_output(self, mocked_print):
        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            metric="nrf",
            clusters=None,
            n_pcs=None,
            top_loadings=5,
            plot=None,
            json=True,
        )
        svc = SpectralDiscordance(args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        assert "scores" in payload
        assert "variance_explained" in payload
        assert "cluster_assignments" in payload
        assert "top_loadings" in payload
        assert "eigengap_values" in payload
        assert "n_clusters" in payload
        assert "metric" in payload

    @patch("builtins.print")
    def test_no_species_tree(self, mocked_print, no_species_tree_args):
        svc = SpectralDiscordance(no_species_tree_args)
        svc.run()
        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Spectral Discordance" in all_output

    @patch("builtins.print")
    def test_wrf_metric(self, mocked_print, wrf_args):
        svc = SpectralDiscordance(wrf_args)
        svc.run()
        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "wrf" in all_output.lower()

    @patch("builtins.print")
    def test_json_loadings_have_species_tree_flag(self, mocked_print):
        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            metric="nrf",
            clusters=None,
            n_pcs=None,
            top_loadings=5,
            plot=None,
            json=True,
        )
        svc = SpectralDiscordance(args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        pc1_loadings = payload["top_loadings"]["PC1"]
        assert len(pc1_loadings) > 0
        assert "in_species_tree" in pc1_loadings[0]
        assert "bipartition" in pc1_loadings[0]
        assert "loading" in pc1_loadings[0]


class TestPlot:
    @patch("builtins.print")
    def test_plots_created(self, mocked_print):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = os.path.join(tmpdir, "test")
            args = Namespace(
                tree=TREE_SIMPLE,
                gene_trees=GENE_TREES,
                metric="nrf",
                clusters=None,
                n_pcs=None,
                top_loadings=5,
                plot=prefix,
                json=False,
            )
            svc = SpectralDiscordance(args)
            svc.run()

            scatter = prefix + "_scatter.png"
            eigengap = prefix + "_eigengap.png"
            assert os.path.exists(scatter), f"Missing {scatter}"
            assert os.path.exists(eigengap), f"Missing {eigengap}"
            assert os.path.getsize(scatter) > 0
            assert os.path.getsize(eigengap) > 0


class TestRValidation:
    """Validate mathematical properties that must match R's prcomp().

    Properties validated:
    - Euclidean distances between rows of the centered bipartition matrix
      equal sqrt(RF) for the binary (nrf) case
    - Variance explained sums to 1 and is in descending order
    - Reconstruction: scores @ loadings ≈ X_centered
    - Orthogonality: loadings are orthonormal rows
    """

    def _setup(self):
        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            metric="nrf",
            clusters=None,
            n_pcs=None,
            top_loadings=5,
            plot=None,
            json=False,
        )
        svc = SpectralDiscordance(args)
        trees = svc._parse_gene_trees(GENE_TREES)
        species_tree = svc.read_tree_file()
        shared = svc._get_shared_taxa(trees, species_tree)
        X, bip_index, sp_flags = svc._build_bipartition_matrix(
            trees, shared, metric="nrf", species_tree=species_tree
        )
        scores, ve, loadings = svc._run_pca(X)
        return svc, X, scores, ve, loadings, bip_index, trees, shared

    def test_euclidean_equals_sqrt_rf(self):
        """Euclidean distance in binary bipartition space = sqrt(RF)."""
        svc, X, _, _, _, bip_index, _, _ = self._setup()
        from scipy.spatial.distance import pdist, squareform
        euclidean = squareform(pdist(X, metric="euclidean"))
        manhattan = squareform(pdist(X, metric="cityblock"))  # L1 = RF
        # Euclidean = sqrt(L1) for binary data
        for i in range(X.shape[0]):
            for j in range(i + 1, X.shape[0]):
                expected = np.sqrt(manhattan[i, j])
                assert euclidean[i, j] == pytest.approx(expected, abs=1e-10)

    def test_reconstruction(self):
        """scores @ loadings ≈ X_centered (full rank)."""
        _, X, scores, _, loadings, _, _, _ = self._setup()
        X_centered = X - X.mean(axis=0)
        reconstructed = scores @ loadings
        assert np.allclose(reconstructed, X_centered, atol=1e-10)

    def test_loadings_orthonormal(self):
        """Loading vectors should be orthonormal."""
        _, _, _, _, loadings, _, _, _ = self._setup()
        # Vt @ V = I (rows of loadings are orthonormal)
        product = loadings @ loadings.T
        assert np.allclose(product, np.eye(loadings.shape[0]), atol=1e-10)

    def test_scores_uncorrelated(self):
        """PC scores should be uncorrelated (covariance is diagonal)."""
        _, _, scores, _, _, _, _, _ = self._setup()
        cov = np.cov(scores.T)
        # Off-diagonal should be ~0
        mask = ~np.eye(cov.shape[0], dtype=bool)
        assert np.allclose(cov[mask], 0, atol=1e-10)
