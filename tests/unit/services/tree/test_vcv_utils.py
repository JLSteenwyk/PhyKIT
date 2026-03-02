import pytest
import numpy as np
from io import StringIO
from pathlib import Path

from Bio import Phylo

from phykit.services.tree.vcv_utils import (
    build_vcv_matrix,
    parse_gene_trees,
    build_discordance_vcv,
    _nearest_psd,
)
from phykit.errors import PhykitUserError


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
GENE_TREES_FILE = str(SAMPLE_FILES / "gene_trees_simple.nwk")


def _read_tree(path):
    return Phylo.read(path, "newick")


def _make_tree(newick):
    return Phylo.read(StringIO(newick), "newick")


class TestBuildVcvMatrix:
    def test_symmetric(self):
        tree = _read_tree(TREE_SIMPLE)
        tips = sorted(t.name for t in tree.get_terminals())
        vcv = build_vcv_matrix(tree, tips)
        np.testing.assert_array_almost_equal(vcv, vcv.T)

    def test_diagonal_is_root_to_tip(self):
        tree = _read_tree(TREE_SIMPLE)
        tips = sorted(t.name for t in tree.get_terminals())
        vcv = build_vcv_matrix(tree, tips)
        for i, name in enumerate(tips):
            expected = tree.distance(tree.root, name)
            assert vcv[i, i] == pytest.approx(expected, rel=1e-6)

    def test_correct_shape(self):
        tree = _read_tree(TREE_SIMPLE)
        tips = sorted(t.name for t in tree.get_terminals())
        vcv = build_vcv_matrix(tree, tips)
        assert vcv.shape == (8, 8)

    def test_positive_semi_definite(self):
        tree = _read_tree(TREE_SIMPLE)
        tips = sorted(t.name for t in tree.get_terminals())
        vcv = build_vcv_matrix(tree, tips)
        eigvals = np.linalg.eigvalsh(vcv)
        assert np.all(eigvals >= -1e-10)

    def test_matches_existing_service_method(self):
        """Shared utility must produce identical output to old inline method."""
        from phykit.services.tree.phylogenetic_signal import PhylogeneticSignal
        from argparse import Namespace

        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=str(SAMPLE_FILES / "tree_simple_traits.tsv"),
        )
        svc = PhylogeneticSignal(args)
        tree = svc.read_tree_file()
        tips = sorted(svc.get_tip_names_from_tree(tree))

        vcv_service = svc._build_vcv_matrix(tree, tips)
        vcv_utility = build_vcv_matrix(tree, tips)
        np.testing.assert_array_almost_equal(vcv_service, vcv_utility)

    def test_simple_three_taxon_tree(self):
        tree = _make_tree("(A:1.0,(B:0.5,C:0.5):0.5);")
        tips = ["A", "B", "C"]
        vcv = build_vcv_matrix(tree, tips)
        # A: root-to-tip = 1.0
        assert vcv[0, 0] == pytest.approx(1.0, rel=1e-6)
        # B: root-to-tip = 1.0
        assert vcv[1, 1] == pytest.approx(1.0, rel=1e-6)
        # C: root-to-tip = 1.0
        assert vcv[2, 2] == pytest.approx(1.0, rel=1e-6)
        # B-C shared path = 0.5
        assert vcv[1, 2] == pytest.approx(0.5, rel=1e-6)
        # A-B shared path = 0 (diverge at root)
        assert vcv[0, 1] == pytest.approx(0.0, abs=1e-6)


class TestParseGeneTrees:
    def test_parse_multi_newick(self):
        trees = parse_gene_trees(GENE_TREES_FILE)
        assert len(trees) == 10

    def test_each_tree_has_tips(self):
        trees = parse_gene_trees(GENE_TREES_FILE)
        for t in trees:
            tips = [tip.name for tip in t.get_terminals()]
            assert len(tips) == 8

    def test_file_not_found(self):
        with pytest.raises(PhykitUserError):
            parse_gene_trees("/nonexistent/path.nwk")

    def test_empty_file(self, tmp_path):
        empty = tmp_path / "empty.nwk"
        empty.write_text("")
        with pytest.raises(PhykitUserError):
            parse_gene_trees(str(empty))

    def test_comments_only(self, tmp_path):
        f = tmp_path / "comments.nwk"
        f.write_text("# comment 1\n# comment 2\n")
        with pytest.raises(PhykitUserError):
            parse_gene_trees(str(f))

    def test_single_tree(self, tmp_path):
        f = tmp_path / "single.nwk"
        f.write_text("(A:1.0,B:2.0,C:3.0);\n")
        trees = parse_gene_trees(str(f))
        assert len(trees) == 1

    def test_comments_ignored(self, tmp_path):
        f = tmp_path / "with_comments.nwk"
        f.write_text("# header\n(A:1.0,B:2.0,C:3.0);\n# footer\n")
        trees = parse_gene_trees(str(f))
        assert len(trees) == 1


class TestBuildDiscordanceVcv:
    def test_identical_gene_trees_match_single_tree(self):
        """When all gene trees are identical to species tree, averaged VCV
        should match the species tree VCV."""
        tree = _read_tree(TREE_SIMPLE)
        tips = sorted(t.name for t in tree.get_terminals())

        # Create 3 copies of the species tree as gene trees
        gene_trees = [_read_tree(TREE_SIMPLE) for _ in range(3)]

        vcv_species = build_vcv_matrix(tree, tips)
        vcv_disc, meta = build_discordance_vcv(tree, gene_trees, tips)

        np.testing.assert_array_almost_equal(vcv_disc, vcv_species, decimal=5)
        assert meta["n_gene_trees"] == 3
        assert meta["n_shared_taxa"] == 8
        assert meta["psd_corrected"] is False

    def test_discordant_trees_differ_from_species_tree(self):
        """When gene trees have different topologies, the averaged VCV
        should differ from the species tree VCV."""
        tree = _read_tree(TREE_SIMPLE)
        gene_trees = parse_gene_trees(GENE_TREES_FILE)
        tips = sorted(t.name for t in tree.get_terminals())

        vcv_species = build_vcv_matrix(tree, tips)
        vcv_disc, meta = build_discordance_vcv(tree, gene_trees, tips)

        # They should not be identical (gene trees 7-9 have discordant topologies)
        assert not np.allclose(vcv_disc, vcv_species, atol=0.1)
        assert meta["n_gene_trees"] == 10

    def test_result_is_psd(self):
        tree = _read_tree(TREE_SIMPLE)
        gene_trees = parse_gene_trees(GENE_TREES_FILE)
        tips = sorted(t.name for t in tree.get_terminals())

        vcv_disc, _ = build_discordance_vcv(tree, gene_trees, tips)

        eigvals = np.linalg.eigvalsh(vcv_disc)
        assert np.all(eigvals >= -1e-10)

    def test_result_is_symmetric(self):
        tree = _read_tree(TREE_SIMPLE)
        gene_trees = parse_gene_trees(GENE_TREES_FILE)
        tips = sorted(t.name for t in tree.get_terminals())

        vcv_disc, _ = build_discordance_vcv(tree, gene_trees, tips)
        np.testing.assert_array_almost_equal(vcv_disc, vcv_disc.T)

    def test_taxa_intersection(self):
        """Gene trees with subset of taxa should auto-prune."""
        tree = _read_tree(TREE_SIMPLE)
        # Create a gene tree with only 5 taxa
        gt = _make_tree("(raccoon:10,(bear:5,(cat:3,dog:3):2):5);")
        gene_trees = [gt]
        all_tips = sorted(t.name for t in tree.get_terminals())

        vcv_disc, meta = build_discordance_vcv(tree, gene_trees, all_tips)
        assert meta["n_shared_taxa"] == 4  # raccoon, bear, cat, dog
        assert set(meta["shared_taxa"]) == {"bear", "cat", "dog", "raccoon"}

    def test_too_few_shared_taxa(self):
        tree = _read_tree(TREE_SIMPLE)
        gt = _make_tree("(taxon_x:1.0,taxon_y:2.0);")
        with pytest.raises(PhykitUserError):
            build_discordance_vcv(tree, [gt], ["taxon_x", "taxon_y"])

    def test_gene_tree_without_branch_lengths(self):
        tree = _read_tree(TREE_SIMPLE)
        gt = _make_tree("(raccoon,(bear,(cat,dog)));")
        tips = sorted(t.name for t in tree.get_terminals())
        with pytest.raises(PhykitUserError):
            build_discordance_vcv(tree, [gt], tips)

    def test_single_gene_tree(self):
        """Single gene tree should produce that tree's VCV."""
        gt = _make_tree("(A:1.0,(B:0.5,C:0.5):0.5);")
        species = _make_tree("(A:2.0,(B:1.0,C:1.0):1.0);")

        gt_tips = sorted(t.name for t in gt.get_terminals())

        vcv_single, meta = build_discordance_vcv(species, [gt], gt_tips)
        vcv_gt = build_vcv_matrix(gt, gt_tips)
        np.testing.assert_array_almost_equal(vcv_single, vcv_gt)
        assert meta["n_gene_trees"] == 1

    def test_metadata_contents(self):
        tree = _read_tree(TREE_SIMPLE)
        gene_trees = parse_gene_trees(GENE_TREES_FILE)
        tips = sorted(t.name for t in tree.get_terminals())

        _, meta = build_discordance_vcv(tree, gene_trees, tips)
        assert "n_gene_trees" in meta
        assert "n_shared_taxa" in meta
        assert "shared_taxa" in meta
        assert "psd_corrected" in meta
        assert "min_eigenvalue_pre_correction" in meta


class TestNearestPsd:
    def test_identity_unchanged(self):
        I = np.eye(4)
        result, corrected, min_eval = _nearest_psd(I)
        np.testing.assert_array_almost_equal(result, I)
        assert corrected is False
        assert min_eval == pytest.approx(1.0)

    def test_already_psd_unchanged(self):
        tree = _read_tree(TREE_SIMPLE)
        tips = sorted(t.name for t in tree.get_terminals())
        vcv = build_vcv_matrix(tree, tips)
        result, corrected, _ = _nearest_psd(vcv)
        np.testing.assert_array_almost_equal(result, vcv)
        assert corrected is False

    def test_negative_eigenvalue_corrected(self):
        # Create a matrix with a negative eigenvalue
        A = np.array([
            [2.0, 1.5, 0.5],
            [1.5, 2.0, 1.5],
            [0.5, 1.5, 2.0],
        ])
        # Force a negative eigenvalue by subtracting from the smallest eigval direction
        eigvals, eigvecs = np.linalg.eigh(A)
        eigvals[0] = -0.5
        bad = eigvecs @ np.diag(eigvals) @ eigvecs.T

        result, corrected, min_eval = _nearest_psd(bad)
        assert corrected is True
        assert min_eval < 0

        # Result should be PSD
        result_eigvals = np.linalg.eigvalsh(result)
        assert np.all(result_eigvals >= -1e-10)

    def test_result_is_symmetric(self):
        A = np.array([
            [2.0, 1.5, 0.5],
            [1.5, 2.0, 1.5],
            [0.5, 1.5, 2.0],
        ])
        eigvals, eigvecs = np.linalg.eigh(A)
        eigvals[0] = -0.5
        bad = eigvecs @ np.diag(eigvals) @ eigvecs.T

        result, _, _ = _nearest_psd(bad)
        np.testing.assert_array_almost_equal(result, result.T)


class TestBackwardCompatibility:
    """Each service should produce identical output when no gene trees are provided."""

    def test_signal_no_gene_trees(self):
        from phykit.services.tree.phylogenetic_signal import PhylogeneticSignal
        from argparse import Namespace

        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=str(SAMPLE_FILES / "tree_simple_traits.tsv"),
            method="blombergs_k",
            permutations=100,
            json=True,
        )
        svc = PhylogeneticSignal(args)
        assert svc.gene_trees_path is None
        # Should run without error
        import io, sys
        old_stdout = sys.stdout
        sys.stdout = io.StringIO()
        try:
            svc.run()
        finally:
            sys.stdout = old_stdout

    def test_regression_no_gene_trees(self):
        from phykit.services.tree.phylogenetic_regression import PhylogeneticRegression
        from argparse import Namespace

        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=str(SAMPLE_FILES / "tree_simple_multi_traits.tsv"),
            response="brain_size",
            predictors=["body_mass"],
            method="BM",
            json=True,
        )
        svc = PhylogeneticRegression(args)
        assert svc.gene_trees_path is None
        import io, sys
        old_stdout = sys.stdout
        sys.stdout = io.StringIO()
        try:
            svc.run()
        finally:
            sys.stdout = old_stdout

    def test_fit_continuous_no_gene_trees(self):
        from phykit.services.tree.fit_continuous import FitContinuous
        from argparse import Namespace

        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=str(SAMPLE_FILES / "tree_simple_traits.tsv"),
            models="BM",
            json=True,
        )
        svc = FitContinuous(args)
        assert svc.gene_trees_path is None
        import io, sys
        old_stdout = sys.stdout
        sys.stdout = io.StringIO()
        try:
            svc.run()
        finally:
            sys.stdout = old_stdout
