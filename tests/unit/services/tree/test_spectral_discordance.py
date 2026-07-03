import json
import os
import subprocess
import sys
import tempfile

import pytest
import numpy as np
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch

from phykit.services.tree.spectral_discordance import SpectralDiscordance
import phykit.services.tree.spectral_discordance as spectral_discordance_module
from phykit.errors import PhykitUserError

here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
GENE_TREES = str(SAMPLE_FILES / "gene_trees_simple.nwk")


def test_module_import_does_not_import_numpy_biophylo_scipy_or_matplotlib():
    code = """
import sys
import phykit.services.tree.spectral_discordance as module

assert hasattr(module.np, "__getattr__")
assert callable(module.print_json)
assert "typing" not in sys.modules
assert "json" not in sys.modules
assert "pickle" not in sys.modules
assert "numpy" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "scipy" not in sys.modules
assert "matplotlib" not in sys.modules
assert "matplotlib.pyplot" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_lazy_numpy_proxy_caches_resolved_attributes():
    lazy_np = spectral_discordance_module._LazyNumpy()

    first_argmin = lazy_np.argmin
    second_argmin = lazy_np.argmin

    assert first_argmin is second_argmin
    assert lazy_np.__dict__["argmin"] is first_argmin


def test_shared_gene_tree_taxa_does_not_slice_gene_trees():
    class NoSliceList(list):
        def __getitem__(self, key):
            if isinstance(key, slice):
                raise AssertionError("gene-tree shared taxa scan should not slice")
            return super().__getitem__(key)

    trees = NoSliceList(["t0", "t1", "t2"])
    tip_sets = {
        "t0": ("A", "B", "C", "D"),
        "t1": ("A", "B", "C"),
        "t2": ("B", "C", "E"),
    }

    assert spectral_discordance_module._shared_gene_tree_taxa(
        trees,
        tip_sets.__getitem__,
    ) == {"B", "C"}


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

    def test_parse_gene_tree_path_list_avoids_per_row_path_objects(
        self, tmp_path, default_args, mocker
    ):
        tree_a = tmp_path / "a.nwk"
        tree_b = tmp_path / "b.nwk"
        tree_a.write_text("((A,B),(C,D));\n")
        tree_b.write_text("((A,C),(B,D));\n")
        tree_list = tmp_path / "tree-list.txt"
        tree_list.write_text(f"{tree_a}\n{tree_b}\n")
        path_calls = 0
        original_path = spectral_discordance_module.Path

        def counting_path(*args, **kwargs):
            nonlocal path_calls
            path_calls += 1
            return original_path(*args, **kwargs)

        mocker.patch.object(
            spectral_discordance_module,
            "Path",
            side_effect=counting_path,
        )
        svc = SpectralDiscordance(default_args)

        trees = svc._parse_gene_trees(str(tree_list))

        assert len(trees) == 2
        assert path_calls == 1

    def test_shared_taxa(self, default_args):
        svc = SpectralDiscordance(default_args)
        trees = svc._parse_gene_trees(GENE_TREES)
        shared = svc._get_shared_taxa(trees, species_tree=None)
        assert len(shared) == 8

    def test_shared_taxa_uses_fast_tip_name_helper(self, default_args, mocker):
        svc = SpectralDiscordance(default_args)
        trees = svc._parse_gene_trees(GENE_TREES)
        species_tree = svc.read_tree_file()
        spy = mocker.spy(svc, "get_tip_names_from_tree")

        shared = svc._get_shared_taxa(trees, species_tree=species_tree)

        assert len(shared) == 8
        assert spy.call_count == len(trees) + 1


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

    def test_canonical_split_equal_size_uses_sorted_tiebreak(self, default_args):
        svc = SpectralDiscordance(default_args)
        all_taxa = frozenset({"A", "B", "C", "D"})
        result = svc._canonical_split(frozenset({"C", "D"}), all_taxa)
        assert result == frozenset({"A", "B"})

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

    def test_extract_splits_uses_direct_postorder(self, default_args):
        from Bio import Phylo
        from io import StringIO

        svc = SpectralDiscordance(default_args)
        tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1,E:1):1);"), "newick"
        )
        all_taxa = frozenset(["A", "B", "C", "D"])

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("standard tree split extraction should use direct postorder")

        tree.find_clades = fail_find_clades

        assert svc._extract_splits(tree, all_taxa) == {frozenset({"A", "B"})}

    def test_extract_splits_with_lengths_uses_direct_postorder(self, default_args):
        from Bio import Phylo
        from io import StringIO

        svc = SpectralDiscordance(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):2,(C:1,D:1):3);"), "newick")
        all_taxa = frozenset(["A", "B", "C", "D"])

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("standard tree split extraction should use direct postorder")

        tree.find_clades = fail_find_clades

        assert svc._extract_splits_with_lengths(tree, all_taxa) == {
            frozenset({"A", "B"}): 3.0
        }

    def test_extract_splits_handles_trifurcating_root_children(self, default_args):
        from Bio import Phylo
        from io import StringIO

        svc = SpectralDiscordance(default_args)
        tree = Phylo.read(
            StringIO("((A:1,B:1):2,(C:1,D:1):3,(E:1,F:1):4);"),
            "newick",
        )
        all_taxa = frozenset(["A", "B", "C", "D", "E", "F"])

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("standard tree split extraction should use direct postorder")

        tree.find_clades = fail_find_clades

        expected_splits = {
            frozenset({"A", "B"}),
            frozenset({"C", "D"}),
            frozenset({"E", "F"}),
        }
        assert svc._extract_splits(tree, all_taxa) == expected_splits
        assert svc._extract_splits_with_lengths(tree, all_taxa) == {
            frozenset({"A", "B"}): 2.0,
            frozenset({"C", "D"}): 3.0,
            frozenset({"E", "F"}): 4.0,
        }

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

    def test_bipartition_matrix_skips_copy_when_no_pruning_needed(self, default_args):
        from Bio import Phylo
        from io import StringIO

        svc = SpectralDiscordance(default_args)
        trees = [
            Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick"),
            Phylo.read(StringIO("((A:1,C:1):1,(B:1,D:1):1);"), "newick"),
        ]
        shared = {"A", "B", "C", "D"}

        with patch(
            "phykit.services.tree.spectral_discordance.pickle.dumps",
            side_effect=AssertionError("unexpected tree copy"),
        ):
            X, bip_index, sp_flags = svc._build_bipartition_matrix(
                trees,
                shared,
                metric="nrf",
                species_tree=trees[0],
            )

        assert X.shape[0] == 2
        assert X.shape[1] == len(bip_index)
        assert len(sp_flags) == len(bip_index)

    def test_copy_prune_if_needed_no_prune_uses_direct_traversal(
        self, default_args, monkeypatch
    ):
        from Bio import Phylo
        from Bio.Phylo.BaseTree import TreeMixin
        from io import StringIO

        svc = SpectralDiscordance(default_args)
        tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1,E:1):1);"), "newick"
        )

        def fail_traversal(*_args, **_kwargs):
            raise AssertionError("generic traversal should not be used")

        def fail_copy(*_args, **_kwargs):
            raise AssertionError("tree should not be copied")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)
        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)
        monkeypatch.setattr(
            "phykit.services.tree.spectral_discordance.pickle.dumps",
            fail_copy,
        )

        pruned = svc._copy_prune_if_needed(tree, {"A", "B", "C", "D", "E"})

        assert pruned is tree

    def test_copy_prune_if_needed_does_not_mutate_original(self, default_args, monkeypatch):
        from Bio import Phylo
        from Bio.Phylo.BaseTree import TreeMixin
        from io import StringIO

        svc = SpectralDiscordance(default_args)
        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1,G:1);"),
            "newick",
        )

        def fail_prune(*_args, **_kwargs):
            raise AssertionError("per-target prune should not be used")

        monkeypatch.setattr(TreeMixin, "prune", fail_prune)

        pruned = svc._copy_prune_if_needed(tree, {"A", "B", "D", "F"})

        assert {tip.name for tip in tree.get_terminals()} == {
            "A", "B", "C", "D", "E", "F", "G",
        }
        assert {tip.name for tip in pruned.get_terminals()} == {"A", "B", "D", "F"}


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

    def test_singular_value_total_variance_uses_dot_product(self, monkeypatch):
        singular_values = np.array([1.0, 2.0, 4.0, 8.0])
        expected = float(np.sum(singular_values ** 2))

        monkeypatch.setattr(
            spectral_discordance_module.np,
            "sum",
            lambda *args, **kwargs: pytest.fail(
                "singular value variance should use np.dot"
            ),
        )

        assert spectral_discordance_module._singular_value_total_variance(
            singular_values
        ) == pytest.approx(expected)

    def test_row_l2_norms_matches_linalg_norm(self):
        matrix = np.array(
            [
                [3.0, 4.0, 0.0],
                [1.0, 2.0, 2.0],
                [0.0, 0.0, 0.0],
            ]
        )

        observed = spectral_discordance_module._row_l2_norms(matrix)
        expected = np.linalg.norm(matrix, axis=1, keepdims=True)

        np.testing.assert_allclose(observed, expected)

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

    def test_top_loading_indices_uses_partial_selection_without_full_argsort(
        self, monkeypatch
    ):
        loadings = np.array([0.1, -0.8, 0.3, 0.6, -0.2])

        def fail_argsort(*_args, **_kwargs):
            raise AssertionError("unique top loadings should not use full argsort")

        monkeypatch.setattr(
            "phykit.services.tree.spectral_discordance.np.argsort",
            fail_argsort,
        )

        top_idx = SpectralDiscordance._top_loading_indices(loadings, 2)

        np.testing.assert_array_equal(top_idx, np.array([1, 3]))

    def test_top_loading_indices_boundary_ties_match_full_sort(self):
        loadings = np.array([0.1, -0.9, 0.9, 0.7])
        expected = np.argsort(np.abs(loadings))[::-1][:1]

        top_idx = SpectralDiscordance._top_loading_indices(loadings, 1)

        np.testing.assert_array_equal(top_idx, expected)


class TestSpectralClustering:
    @staticmethod
    def _reference_kmeans(X: np.ndarray, K: int, max_iter: int = 300) -> np.ndarray:
        G = X.shape[0]
        rng = np.random.RandomState(42)

        centers = [X[rng.randint(G)]]
        for _ in range(1, K):
            dists = np.array([
                min(np.sum((x - c) ** 2) for c in centers)
                for x in X
            ])
            probs = dists / dists.sum() if dists.sum() > 0 else np.ones(G) / G
            centers.append(X[rng.choice(G, p=probs)])
        centers = np.array(centers)

        for _ in range(max_iter):
            dists_to_centers = np.array([
                np.sum((X - c) ** 2, axis=1) for c in centers
            ]).T
            labels = np.argmin(dists_to_centers, axis=1)

            new_centers = np.array([
                X[labels == k].mean(axis=0) if np.any(labels == k) else centers[k]
                for k in range(K)
            ])
            if np.allclose(new_centers, centers):
                break
            centers = new_centers

        return labels

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

    def test_squared_distances_to_centers_matches_reference(self):
        X = np.array([[0.0, 1.0], [2.0, 3.0], [4.0, 5.0]])
        centers = np.array([[0.0, 0.0], [1.0, 1.0]])

        observed = spectral_discordance_module._squared_distances_to_centers(
            X,
            centers,
        )
        expected = np.array([
            [1.0, 1.0],
            [13.0, 5.0],
            [41.0, 25.0],
        ])

        np.testing.assert_allclose(observed, expected)

    def test_squared_distances_to_center_matches_reference(self):
        X = np.array([[0.0, 1.0], [2.0, 3.0], [4.0, 5.0]])
        center = np.array([1.0, 1.0])

        observed = spectral_discordance_module._squared_distances_to_center(
            X,
            center,
        )

        np.testing.assert_allclose(observed, np.array([1.0, 5.0, 25.0]))

    def test_update_kmeans_centers_preserves_empty_clusters(self):
        X = np.array([[0.0, 0.0], [2.0, 2.0], [4.0, 4.0]])
        labels = np.array([0, 0, 2])
        centers = np.array([[10.0, 10.0], [20.0, 20.0], [30.0, 30.0]])

        updated = spectral_discordance_module._update_kmeans_centers(
            X,
            labels,
            centers,
        )

        np.testing.assert_allclose(
            updated,
            np.array([[1.0, 1.0], [20.0, 20.0], [4.0, 4.0]]),
        )

    def test_kmeans_matches_reference_labels(self):
        rng = np.random.default_rng(20260628)
        X = rng.normal(size=(200, 5))

        expected = self._reference_kmeans(X, 4, max_iter=20)
        observed = SpectralDiscordance._kmeans(X, 4, max_iter=20)

        np.testing.assert_array_equal(observed, expected)

    def test_kmeans_initialization_uses_einsum_distances(self, monkeypatch):
        rng = np.random.default_rng(20260701)
        X = rng.normal(size=(200, 5))

        def fail_sum(*_args, **_kwargs):
            raise AssertionError("kmeans initialization should avoid np.sum")

        monkeypatch.setattr(spectral_discordance_module.np, "sum", fail_sum)

        labels = SpectralDiscordance._kmeans(X, 4, max_iter=20)

        assert labels.shape == (200,)

    def test_override_k(self, default_args):
        labels, K, _ = self._get_clusters(default_args, n_clusters=3)
        assert K == 3
        assert len(set(labels)) <= 3

    def test_spectral_cluster_matches_dense_laplacian_reference(
        self, default_args, monkeypatch
    ):
        import phykit.services.tree.spectral_discordance as sd
        from scipy.spatial.distance import pdist, squareform

        svc = SpectralDiscordance(default_args)
        rng = np.random.default_rng(20260623)
        X_centered = rng.normal(size=(18, 5))
        dists = squareform(pdist(X_centered, metric="euclidean"))
        upper_tri = dists[np.triu_indices(X_centered.shape[0], k=1)]
        sigma = float(np.median(upper_tri))
        W = np.exp(-(dists ** 2) / (2 * sigma ** 2))
        np.fill_diagonal(W, W.diagonal() + 1e-10)
        d_inv_sqrt = 1.0 / np.sqrt(np.sum(W, axis=1))
        D_inv_sqrt = np.diag(d_inv_sqrt)
        dense_laplacian = np.eye(X_centered.shape[0]) - D_inv_sqrt @ W @ D_inv_sqrt
        eigenvalues = np.linalg.eigvalsh(dense_laplacian)
        max_k = min(20, X_centered.shape[0] // 2)
        expected_eigengaps = np.diff(eigenvalues[:max_k + 1])
        expected_k = int(np.argmax(expected_eigengaps[1:max_k]) + 2)
        expected_k = max(expected_k, 2)

        original_diag = sd.np.diag

        def diag_without_vector_to_matrix(a, *args, **kwargs):
            arr = np.asarray(a)
            if arr.ndim == 1:
                raise AssertionError("spectral clustering should not allocate D_inv")
            return original_diag(a, *args, **kwargs)

        monkeypatch.setattr(sd.np, "diag", diag_without_vector_to_matrix)
        monkeypatch.setattr(
            sd.np,
            "sum",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("spectral clustering should use ndarray row sums")
            ),
        )
        labels, K, eigengaps = svc._spectral_cluster(X_centered)

        assert len(labels) == X_centered.shape[0]
        assert K == expected_k
        np.testing.assert_allclose(eigengaps, expected_eigengaps)

    def test_spectral_cluster_reuses_condensed_distances(
        self, default_args, monkeypatch
    ):
        svc = SpectralDiscordance(default_args)
        X_centered = np.array(
            [
                [0.0, 0.0],
                [0.1, 0.0],
                [0.0, 0.1],
                [4.0, 4.0],
                [4.1, 4.0],
                [4.0, 4.1],
            ]
        )

        def fail_triu_indices(*_args, **_kwargs):
            raise AssertionError(
                "spectral clustering should reuse pdist's condensed distances"
            )

        monkeypatch.setattr(
            spectral_discordance_module.np,
            "triu_indices",
            fail_triu_indices,
        )

        labels, K, eigengaps = svc._spectral_cluster(X_centered, n_clusters=2)

        assert len(labels) == X_centered.shape[0]
        assert K == 2
        assert eigengaps.size > 0

    def test_spectral_cluster_row_normalization_avoids_linalg_norm(
        self, default_args, monkeypatch
    ):
        svc = SpectralDiscordance(default_args)
        X_centered = np.array(
            [
                [0.0, 0.0],
                [0.1, 0.0],
                [0.0, 0.1],
                [4.0, 4.0],
                [4.1, 4.0],
                [4.0, 4.1],
            ]
        )

        def fail_norm(*_args, **_kwargs):
            raise AssertionError(
                "spectral clustering should use the row L2 helper"
            )

        monkeypatch.setattr(
            spectral_discordance_module.np.linalg,
            "norm",
            fail_norm,
        )

        labels, K, eigengaps = svc._spectral_cluster(X_centered, n_clusters=2)

        assert len(labels) == X_centered.shape[0]
        assert K == 2
        assert eigengaps.size > 0

    def test_concordant_trees_cluster_together(self, default_args):
        labels, K, _ = self._get_clusters(default_args)
        assert K >= 2
        concordant_labels = [labels[i] for i in range(6)]
        from collections import Counter
        most_common = Counter(concordant_labels).most_common(1)[0][1]
        assert most_common >= 4


class TestRun:
    def test_run_reuses_unmodified_species_tree(self, default_args):
        svc = SpectralDiscordance(default_args)
        species_tree = object()
        gene_trees = [object()] * 5
        shared = {"A", "B", "C", "D"}
        X = np.zeros((5, 3))
        bip_index = [frozenset({"A", "B"})]
        sp_flags = {}
        scores = np.zeros((5, 2))
        var_explained = np.array([1.0, 0.0])
        loadings = np.zeros((2, 3))
        labels = np.zeros(5, dtype=int)
        eigengaps = np.array([0.0, 0.0])
        top_loadings = {"PC1": []}

        with patch.object(
            svc,
            "read_tree_file",
            side_effect=AssertionError("run should not copy the cached species tree"),
        ), patch.object(
            svc, "read_tree_file_unmodified", return_value=species_tree
        ) as read_unmodified, patch.object(
            svc, "_parse_gene_trees", return_value=gene_trees
        ), patch.object(
            svc, "_get_shared_taxa", return_value=shared
        ) as get_shared, patch.object(
            svc, "_build_bipartition_matrix",
            return_value=(X, bip_index, sp_flags),
        ) as build_matrix, patch.object(
            svc, "_run_pca", return_value=(scores, var_explained, loadings)
        ), patch.object(
            svc, "_spectral_cluster", return_value=(labels, 2, eigengaps)
        ), patch.object(
            svc, "_get_top_loadings", return_value=top_loadings
        ), patch.object(
            svc, "_print_text"
        ):
            svc.run()

        read_unmodified.assert_called_once_with()
        get_shared.assert_called_once_with(gene_trees, species_tree)
        build_matrix.assert_called_once_with(
            gene_trees, shared, metric="nrf", species_tree=species_tree
        )

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
    def test_text_output_batches_pc_score_rows(self, mocked_print, default_args):
        svc = SpectralDiscordance(default_args)
        scores = np.array(
            [
                [0.1, 0.2, 0.3, 0.4, 0.5],
                [1.1, 1.2, 1.3, 1.4, 1.5],
            ]
        )
        var_explained = np.array([0.4, 0.2, 0.1, 0.05, 0.01])
        top_loadings = {
            "PC1": [
                {
                    "loading": 0.25,
                    "bipartition": "A|B",
                    "in_species_tree": True,
                }
            ]
        }
        labels = np.array([0, 1])

        svc._print_text(
            scores,
            var_explained,
            top_loadings,
            labels,
            2,
            np.array([]),
            5,
            2,
            4,
            6,
        )

        mocked_print.assert_called_once()
        text = mocked_print.call_args.args[0]
        assert "Spectral Discordance Decomposition" in text
        assert "gene_tree_1" in text
        assert "gene_tree_2" in text
        assert "    0.5000" in text
        assert text.count("gene_tree_") == 2

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

    def test_format_json_preserves_payload_shape(self, default_args):
        svc = SpectralDiscordance(default_args)
        scores = np.array([
            [0.1, 0.2, 0.3],
            [1.1, 1.2, 1.3],
        ])
        labels = np.array([1, 0])
        top_loadings = {
            "PC1": [
                {
                    "bipartition": "A|B",
                    "loading": 0.75,
                    "in_species_tree": True,
                }
            ]
        }

        payload = svc._format_json(
            scores,
            np.array([0.6, 0.3, 0.1]),
            top_loadings,
            labels,
            2,
            np.array([0.4, 0.2]),
            2,
            {frozenset({"A", "B"}): 0, frozenset({"C", "D"}): 1},
            {frozenset({"A", "B"}): True},
        )

        assert payload == {
            "metric": "nrf",
            "n_gene_trees": 2,
            "n_bipartitions": 2,
            "n_clusters": 2,
            "scores": {
                "gene_tree_1": {"PC1": 0.1, "PC2": 0.2, "cluster": 1},
                "gene_tree_2": {"PC1": 1.1, "PC2": 1.2, "cluster": 0},
            },
            "variance_explained": {"PC1": 0.6, "PC2": 0.3},
            "top_loadings": top_loadings,
            "cluster_assignments": [1, 0],
            "eigengap_values": [0.4, 0.2],
        }

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
    def test_scatter_plot_reuses_cluster_indices_for_score_coordinates(
        self, tmp_path, default_args
    ):
        try:
            import matplotlib  # noqa: F401
        except ImportError:
            pytest.skip("matplotlib not installed")

        class NoBooleanMaskScores(np.ndarray):
            def __getitem__(self, key):
                if (
                    isinstance(key, tuple)
                    and isinstance(key[0], np.ndarray)
                    and key[0].dtype == bool
                ):
                    raise AssertionError("scatter scores should use cluster indices once")
                return super().__getitem__(key)

        scores = np.asarray(
            [
                [0.0, 0.0],
                [1.0, 1.0],
                [2.0, 4.0],
                [3.0, 9.0],
            ]
        ).view(NoBooleanMaskScores)
        labels = np.asarray([0, 1, 0, 2])

        svc = SpectralDiscordance(default_args)
        with patch("builtins.print"):
            svc._plot_scatter(scores, labels, 3, tmp_path / "scatter.png")

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


class TestPolytomyHandling:
    """Polytomous nodes should be skipped in split extraction."""

    def _make_svc(self):
        return SpectralDiscordance.__new__(SpectralDiscordance)

    def test_star_tree_no_splits(self):
        """A 4-way star (A,B,C,D) produces no splits."""
        from Bio import Phylo
        from io import StringIO
        tree = Phylo.read(StringIO("(A:1,B:1,C:1,D:1);"), "newick")
        all_taxa = frozenset(["A", "B", "C", "D"])
        svc = self._make_svc()
        splits = svc._extract_splits(tree, all_taxa)
        assert len(splits) == 0

    def test_trifurcating_root_allowed(self):
        """A trifurcating root (standard unrooted) is not skipped."""
        from Bio import Phylo
        from io import StringIO
        tree = Phylo.read(StringIO("((A:1,B:1):1,C:1,D:1);"), "newick")
        all_taxa = frozenset(["A", "B", "C", "D"])
        svc = self._make_svc()
        splits = svc._extract_splits(tree, all_taxa)
        assert len(splits) >= 1

    def test_internal_polytomy_skipped(self):
        """An internal polytomy node is excluded from split extraction."""
        from Bio import Phylo
        from io import StringIO
        tree = Phylo.read(StringIO("((A:1,B:1,C:1,D:1):1,E:1,F:1);"), "newick")
        all_taxa = frozenset(["A", "B", "C", "D", "E", "F"])
        svc = self._make_svc()
        splits = svc._extract_splits(tree, all_taxa)
        assert len(splits) == 0

    def test_extract_splits_with_lengths_skips_polytomy(self):
        """_extract_splits_with_lengths also skips polytomous nodes."""
        from Bio import Phylo
        from io import StringIO
        tree = Phylo.read(StringIO("(A:1,B:1,C:1,D:1);"), "newick")
        all_taxa = frozenset(["A", "B", "C", "D"])
        svc = self._make_svc()
        split_lengths = svc._extract_splits_with_lengths(tree, all_taxa)
        assert len(split_lengths) == 0

    def test_resolved_subtree_preserved(self):
        """Resolved subclades still produce splits even with polytomy above."""
        from Bio import Phylo
        from io import StringIO
        tree = Phylo.read(StringIO("(((A:1,B:1):1,C:1,D:1):1,E:1,F:1);"), "newick")
        all_taxa = frozenset(["A", "B", "C", "D", "E", "F"])
        svc = self._make_svc()
        splits = svc._extract_splits(tree, all_taxa)
        assert len(splits) >= 1
