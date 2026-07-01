import json
import os
import subprocess
import sys
import tempfile

import pytest
import numpy as np
from argparse import Namespace
from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin
from io import StringIO
from pathlib import Path
from unittest.mock import patch

from phykit.services.tree.tree_space import TreeSpace
import phykit.services.tree.tree_space as tree_space_module
from phykit.services.tree.base import Tree
from phykit.errors import PhykitUserError


def test_module_import_does_not_import_numpy_or_biophylo():
    code = """
import sys
import phykit.services.tree.tree_space as module
assert hasattr(module.np, "__getattr__")
assert hasattr(module.Phylo, "read")
assert callable(module.print_json)
assert "json" not in sys.modules
assert "pickle" not in sys.modules
assert "typing" not in sys.modules
assert "numpy" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_shared_gene_tree_taxa_does_not_slice_gene_trees():
    class NoSliceList(list):
        def __getitem__(self, key):
            if isinstance(key, slice):
                raise AssertionError("gene-tree shared taxa scan should not slice")
            return super().__getitem__(key)

    trees = NoSliceList(["t0", "t1", "t2"])
    tip_sets = {
        "t0": {"A", "B", "C", "D"},
        "t1": {"A", "B", "C"},
        "t2": {"B", "C", "E"},
    }

    assert tree_space_module._shared_gene_tree_taxa(
        trees,
        tip_sets.__getitem__,
    ) == {"B", "C"}


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


class TestCanonicalSplit:
    def test_equal_size_returns_lexicographically_smaller_side(self):
        all_taxa = frozenset({"A", "B", "C", "D"})

        assert (
            TreeSpace._canonical_split(frozenset({"C", "D"}), all_taxa)
            == frozenset({"A", "B"})
        )
        assert (
            TreeSpace._canonical_split(frozenset({"A", "B"}), all_taxa)
            == frozenset({"A", "B"})
        )

    def test_empty_split_still_canonicalizes_to_empty_set(self):
        assert TreeSpace._canonical_split(frozenset(), frozenset()) == frozenset()


class TestTreeParsing:
    def test_loads_all_gene_trees(self, default_args):
        svc = TreeSpace(default_args)
        trees = svc._parse_trees_from_source(GENE_TREES)
        assert len(trees) == 10

    def test_parse_trees_skips_comments_and_blank_lines(self, tmp_path, default_args):
        tree_file = tmp_path / "trees.nwk"
        tree_file.write_text(
            "# ignored\n\n  ((A,B),(C,D));  \n# also ignored\n((A,C),(B,D));\n"
        )
        svc = TreeSpace(default_args)

        trees = svc._parse_trees_from_source(str(tree_file))

        assert len(trees) == 2
        assert [TreeSpace._tips(tree) for tree in trees] == [
            {"A", "B", "C", "D"},
            {"A", "B", "C", "D"},
        ]

    def test_parse_tree_path_list_avoids_per_row_path_objects(
        self, tmp_path, default_args, mocker
    ):
        tree_a = tmp_path / "a.nwk"
        tree_b = tmp_path / "b.nwk"
        tree_a.write_text("((A,B),(C,D));\n")
        tree_b.write_text("((A,C),(B,D));\n")
        tree_list = tmp_path / "tree-list.txt"
        tree_list.write_text("a.nwk\nb.nwk\n")
        path_calls = 0
        original_path = tree_space_module.Path

        def counting_path(*args, **kwargs):
            nonlocal path_calls
            path_calls += 1
            return original_path(*args, **kwargs)

        mocker.patch.object(tree_space_module, "Path", side_effect=counting_path)
        svc = TreeSpace(default_args)

        trees = svc._parse_trees_from_source(str(tree_list))

        assert len(trees) == 2
        assert path_calls == 1

    def test_parse_tree_path_list_streams_before_consuming_all_rows(
        self, default_args, mocker
    ):
        read_calls = []

        class StreamingHandle:
            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, traceback):
                return False

            def __iter__(self):
                yield "a.nwk\n"
                if not read_calls:
                    raise AssertionError("path rows should parse before full cleanup")
                yield "b.nwk\n"

        class Source:
            parent = Path(".")

            def open(self):
                return StreamingHandle()

        mocker.patch.object(tree_space_module, "Path", return_value=Source())
        mocker.patch.object(tree_space_module, "_path_exists", return_value=True)

        def fake_read(source, _format):
            read_calls.append(source)
            return source

        mocker.patch.object(tree_space_module.Phylo, "read", side_effect=fake_read)
        svc = TreeSpace(default_args)

        assert svc._parse_trees_from_source("tree-list.txt") == ["a.nwk", "b.nwk"]

    def test_file_not_found(self, default_args):
        svc = TreeSpace(default_args)
        with pytest.raises(PhykitUserError):
            svc._parse_trees_from_source("no_such_file.nwk")

    def test_shared_taxa(self, default_args):
        svc = TreeSpace(default_args)
        trees = svc._parse_trees_from_source(GENE_TREES)
        shared = svc._get_shared_taxa(trees)
        assert len(shared) == 8

    def test_shared_taxa_uses_fast_terminal_names_for_parsed_trees(
        self, monkeypatch, default_args
    ):
        svc = TreeSpace(default_args)
        trees = [
            Phylo.read(StringIO("((A,B),(C,D));"), "newick"),
            Phylo.read(StringIO("((A,C),(B,D));"), "newick"),
        ]

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("parsed trees should use fast terminal names")

        for tree in trees:
            monkeypatch.setattr(tree, "get_terminals", fail_get_terminals)

        assert svc._get_shared_taxa(trees) == {"A", "B", "C", "D"}

    def test_prune_to_taxa_uses_standard_tree_batch_prune(
        self, monkeypatch, default_args
    ):
        svc = TreeSpace(default_args)
        tree = Phylo.read(StringIO("(((A,B),(C,D)),((E,F),(G,H)));"), "newick")

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("parsed trees should use fast terminal clades")

        def fail_prune(*args, **kwargs):
            raise AssertionError("standard trees should use batch pruning")

        monkeypatch.setattr(tree, "get_terminals", fail_get_terminals)
        monkeypatch.setattr(tree, "prune", fail_prune)

        pruned = svc._prune_to_taxa(tree, {"A", "C", "E", "G"})

        assert set(Tree.calculate_terminal_names_fast(pruned)) == {
            "A",
            "C",
            "E",
            "G",
        }


class TestDistanceMatrix:
    def test_rf_distance_matrix_symmetric(self, default_args):
        svc = TreeSpace(default_args)
        trees = svc._parse_trees_from_source(GENE_TREES)
        shared = svc._get_shared_taxa(trees)
        dist = svc._build_distance_matrix(trees, shared, "rf")
        np.testing.assert_array_almost_equal(dist, dist.T)

    def test_condensed_distance_vector_matches_upper_triangle_order(self):
        dist = np.array(
            [
                [0.0, 1.0, 2.0, 3.0],
                [1.0, 0.0, 4.0, 5.0],
                [2.0, 4.0, 0.0, 6.0],
                [3.0, 5.0, 6.0, 0.0],
            ],
            dtype=float,
        )

        condensed = tree_space_module._condensed_distance_vector(dist)

        np.testing.assert_array_equal(condensed, np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]))

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

    def test_vectorized_rf_distance_matrix_matches_pairwise_formula(self, default_args):
        svc = TreeSpace(default_args)
        split_sets = [
            {frozenset({"A", "B"}), frozenset({"C", "D"})},
            {frozenset({"A", "B"}), frozenset({"E", "F"})},
            {frozenset({"C", "D"})},
        ]

        dist = svc._build_rf_distance_matrix_from_splits(split_sets, n_taxa=6)

        expected = np.zeros((3, 3))
        for i in range(3):
            for j in range(i + 1, 3):
                expected[i, j] = svc._compute_rf_distance(
                    split_sets[i], split_sets[j], n_taxa=6
                )
                expected[j, i] = expected[i, j]
        np.testing.assert_array_equal(dist, expected)

    def test_vectorized_kf_distance_matrix_matches_pairwise_formula(self, kf_args):
        svc = TreeSpace(kf_args)
        split_lengths = [
            {frozenset({"A", "B"}): 1.0, frozenset({"C", "D"}): 2.0},
            {frozenset({"A", "B"}): 1.5, frozenset({"E", "F"}): 0.5},
            {frozenset({"C", "D"}): 2.25},
        ]

        dist = svc._build_kf_distance_matrix_from_splits(split_lengths)

        expected = np.zeros((3, 3))
        for i in range(3):
            for j in range(i + 1, 3):
                expected[i, j] = svc._compute_kf_distance(
                    split_lengths[i], split_lengths[j]
                )
                expected[j, i] = expected[i, j]
        np.testing.assert_allclose(dist, expected)

    def test_direct_postorder_preserves_left_to_right_order(self, default_args):
        svc = TreeSpace(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):2,(C:1,D:1):3);"), "newick")

        observed = [
            clade.name if clade.name is not None else "internal"
            for clade in svc._postorder_clades_direct(tree)
        ]

        assert observed == [
            "A", "B", "internal", "C", "D", "internal", "internal",
        ]

    def test_extract_splits_uses_direct_postorder_for_standard_tree(
        self, monkeypatch, default_args
    ):
        svc = TreeSpace(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):2,(C:1,D:1):3);"), "newick")

        def fail_find_clades(*args, **kwargs):
            raise AssertionError("standard split extraction should use direct postorder")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        splits = svc._extract_splits(tree, frozenset({"A", "B", "C", "D"}))

        assert splits == {frozenset({"A", "B"})}

    def test_extract_splits_with_lengths_uses_direct_postorder_for_standard_tree(
        self, monkeypatch, kf_args
    ):
        svc = TreeSpace(kf_args)
        tree = Phylo.read(StringIO("((A:1,B:1):2,(C:1,D:1):3);"), "newick")

        def fail_find_clades(*args, **kwargs):
            raise AssertionError("standard split extraction should use direct postorder")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        splits = svc._extract_splits_with_lengths(
            tree,
            frozenset({"A", "B", "C", "D"}),
        )

        assert splits == {frozenset({"A", "B"}): 3.0}

    def test_extract_splits_handles_mixed_child_counts(self, default_args):
        svc = TreeSpace(default_args)
        tree = Phylo.read(
            StringIO("(((A:1):2,((B:3,C:4):5,D:6):7):8,E:9,F:10);"),
            "newick",
        )

        splits = svc._extract_splits(
            tree,
            frozenset({"A", "B", "C", "D", "E", "F"}),
        )

        assert splits == {
            frozenset({"B", "C"}),
            frozenset({"A", "E", "F"}),
            frozenset({"E", "F"}),
        }

    def test_extract_splits_with_lengths_handles_mixed_child_counts(self, kf_args):
        svc = TreeSpace(kf_args)
        tree = Phylo.read(
            StringIO("(((A:1):2,((B:3,C:4):5,D:6):7):8,E:9,F:10);"),
            "newick",
        )

        splits = svc._extract_splits_with_lengths(
            tree,
            frozenset({"A", "B", "C", "D", "E", "F"}),
        )

        assert splits == {
            frozenset({"B", "C"}): 5.0,
            frozenset({"A", "E", "F"}): 7.0,
            frozenset({"E", "F"}): 8.0,
        }

    def test_build_distance_matrix_skips_copy_when_no_pruning_needed(
        self, monkeypatch, default_args
    ):
        svc = TreeSpace(default_args)
        trees = [
            Phylo.read(StringIO("((A,B),(C,D));"), "newick"),
            Phylo.read(StringIO("((A,C),(B,D));"), "newick"),
        ]

        def fail_copy(*args, **kwargs):
            raise AssertionError("tree copy should not be needed")

        monkeypatch.setattr("phykit.services.tree.tree_space.pickle.dumps", fail_copy)

        dist = svc._build_distance_matrix(trees, {"A", "B", "C", "D"}, "rf")

        assert dist.shape == (2, 2)
        np.testing.assert_array_equal(np.diag(dist), 0.0)

    def test_build_distance_matrix_skips_terminal_collection_when_no_pruning_needed(
        self, monkeypatch, default_args
    ):
        svc = TreeSpace(default_args)
        trees = [
            Phylo.read(StringIO("((A,B),(C,D));"), "newick"),
            Phylo.read(StringIO("((A,C),(B,D));"), "newick"),
        ]

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("no-prune setup should use fast terminal names")

        for tree in trees:
            monkeypatch.setattr(tree, "get_terminals", fail_get_terminals)

        dist = svc._build_distance_matrix(trees, {"A", "B", "C", "D"}, "rf")

        assert dist.shape == (2, 2)
        np.testing.assert_array_equal(np.diag(dist), 0.0)

    def test_build_distance_matrix_batch_prunes_copied_terminal_objects(
        self, monkeypatch, default_args
    ):
        svc = TreeSpace(default_args)
        trees = [
            Phylo.read(StringIO("((A,B),(C,D),(E,F));"), "newick"),
            Phylo.read(StringIO("((A,C),(B,D),(E,F));"), "newick"),
        ]

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("copied-prune setup should use direct traversal")

        def fail_prune(*args, **kwargs):
            raise AssertionError("copied-prune setup should use batch pruning")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)
        monkeypatch.setattr(TreeMixin, "prune", fail_prune)

        dist = svc._build_distance_matrix(trees, {"A", "B", "C", "D"}, "rf")

        assert dist.shape == (2, 2)
        np.testing.assert_array_equal(np.diag(dist), 0.0)
        assert set(Tree.calculate_terminal_names_fast(trees[0])) == {
            "A", "B", "C", "D", "E", "F"
        }

    def test_rf_distance_matrix_extracts_shared_splits_without_copying(
        self, monkeypatch, default_args
    ):
        svc = TreeSpace(default_args)
        tree = Phylo.read(StringIO("((A:1,E:1):1,(B:1,(C:1,D:1):1):1);"), "newick")
        shared_taxa = {"A", "B", "C", "D"}
        pruned = Phylo.read(StringIO("((A:1,E:1):1,(B:1,(C:1,D:1):1):1);"), "newick")
        for tip in [tip.name for tip in pruned.get_terminals() if tip.name not in shared_taxa]:
            pruned.prune(tip)
        expected_splits = svc._extract_splits(pruned, frozenset(shared_taxa))

        def fail_copy(*_args, **_kwargs):
            raise AssertionError("RF split extraction should not copy-prune trees")

        monkeypatch.setattr("phykit.services.tree.tree_space.pickle.dumps", fail_copy)

        dist = svc._build_distance_matrix([tree, pruned], shared_taxa, "rf")

        observed_splits = svc._extract_splits(tree, frozenset(shared_taxa))
        assert observed_splits == expected_splits
        assert dist.shape == (2, 2)
        np.testing.assert_array_equal(np.diag(dist), 0.0)


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

    def test_spectral_cluster_uses_condensed_distances(
        self, default_args, monkeypatch
    ):
        svc = TreeSpace(default_args)
        dist = np.array(
            [
                [0.0, 0.1, 0.2, 4.0],
                [0.1, 0.0, 0.3, 4.1],
                [0.2, 0.3, 0.0, 4.2],
                [4.0, 4.1, 4.2, 0.0],
            ]
        )
        original_condensed = tree_space_module._condensed_distance_vector
        calls = []

        def condensed_spy(matrix):
            calls.append(matrix.shape)
            return original_condensed(matrix)

        def fail_triu_indices(*_args, **_kwargs):
            raise AssertionError(
                "spectral clustering should avoid triangular index arrays"
            )

        monkeypatch.setattr(
            tree_space_module,
            "_condensed_distance_vector",
            condensed_spy,
        )
        monkeypatch.setattr(tree_space_module.np, "triu_indices", fail_triu_indices)

        labels, k = svc._spectral_cluster(dist, 2, 42)

        assert calls == [(4, 4)]
        assert len(labels) == 4
        assert k == 2

    def test_auto_k_detection(self, default_args):
        svc = TreeSpace(default_args)
        trees = svc._parse_trees_from_source(GENE_TREES)
        shared = svc._get_shared_taxa(trees)
        dist = svc._build_distance_matrix(trees, shared, "rf")
        labels, k = svc._spectral_cluster(dist, None, 42)
        # k should be a reasonable number
        assert 2 <= k <= len(trees) // 2

    def test_auto_detect_k_matches_dense_laplacian_reference(
        self, default_args, monkeypatch
    ):
        import phykit.services.tree.tree_space as ts

        svc = TreeSpace(default_args)
        rng = np.random.default_rng(20260623)
        points = rng.normal(size=(24, 4))
        distances = np.sqrt(
            np.sum((points[:, None, :] - points[None, :, :]) ** 2, axis=2)
        )
        affinity = np.exp(-(distances ** 2) / (2 * np.median(distances) ** 2))
        np.fill_diagonal(affinity, 1.0)

        d_inv_sqrt = 1.0 / np.sqrt(np.maximum(np.sum(affinity, axis=1), 1e-10))
        D_inv_sqrt = np.diag(d_inv_sqrt)
        dense_laplacian = np.eye(affinity.shape[0]) - (
            D_inv_sqrt @ affinity @ D_inv_sqrt
        )
        eigenvalues = np.sort(np.linalg.eigvalsh(dense_laplacian))
        max_k = min(10, affinity.shape[0] // 2)
        expected_eigengaps = np.diff(eigenvalues[: max_k + 1])
        expected_k = int(np.argmax(expected_eigengaps[1:max_k]) + 2)
        expected_k = max(2, min(expected_k, affinity.shape[0] // 2))

        original_diag = ts.np.diag

        def diag_without_vector_to_matrix(a, *args, **kwargs):
            arr = np.asarray(a)
            if arr.ndim == 1:
                raise AssertionError("auto-K should not allocate D_inv")
            return original_diag(a, *args, **kwargs)

        monkeypatch.setattr(ts.np, "diag", diag_without_vector_to_matrix)
        monkeypatch.setattr(
            ts.np,
            "sum",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("auto-K should use ndarray row sums")
            ),
        )
        assert svc._auto_detect_k(affinity) == expected_k


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

    def test_write_distance_matrix_preserves_csv_format(self, tmp_path):
        output = tmp_path / "dist.csv"
        dist = np.array([
            [0.0, 1.2345678, 2.0],
            [1.2345678, 0.0, 3.5],
            [2.0, 3.5, 0.0],
        ])

        TreeSpace._write_distance_matrix(
            dist,
            ["tree_a", "tree_b", "tree_c"],
            str(output),
        )

        assert output.read_text() == (
            ",tree_a,tree_b,tree_c\n"
            "tree_a,0.000000,1.234568,2.000000\n"
            "tree_b,1.234568,0.000000,3.500000\n"
            "tree_c,2.000000,3.500000,0.000000\n"
        )

    def test_print_text_output_batches_cluster_rows(self, tmp_path, mocker):
        output = tmp_path / "tree_space.png"
        distance_matrix = tmp_path / "dist.csv"
        args = _base_args(
            output=str(output),
            distance_matrix=str(distance_matrix),
            metric="kf",
            method="pcoa",
        )
        svc = TreeSpace(args)
        printed = mocker.patch("builtins.print")

        svc._print_text_output(
            n_gene_trees=7,
            n_taxa=4,
            k=2,
            k_method="user-specified",
            clusters_info=[
                {"id": 1, "n_trees": 3},
                {"id": 2, "n_trees": 4},
            ],
            species_tree_cluster=2,
            plot_mode="heatmap",
        )

        printed.assert_called_once_with(
            "Tree Space Analysis\n"
            "Method: PCOA\n"
            "Metric: KF\n"
            "Gene trees: 7\n"
            "Shared taxa: 4\n"
            "Clusters: 2 (user-specified)\n"
            "  Cluster 1: 3 trees\n"
            "  Cluster 2: 4 trees\n"
            "Species tree: Cluster 2\n"
            "Plot mode: heatmap\n"
            f"Output: {output}\n"
            f"Distance matrix: {distance_matrix}"
        )

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
