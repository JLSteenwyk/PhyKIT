import builtins
import importlib
import subprocess
import sys

import pytest
import numpy as np
from io import StringIO
from pathlib import Path

from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin

import phykit.services.tree.vcv_utils as vcv_utils
from phykit.services.tree.vcv_utils import (
    build_vcv_matrix,
    build_transformed_vcv_matrix,
    parse_gene_trees,
    build_discordance_vcv,
    _build_pruned_subset_vcv_matrix,
    _copy_prune_gene_tree_to_shared_taxa,
    _gene_tree_has_missing_branch_lengths,
    _get_tip_names_and_missing_branch_lengths,
    _nearest_psd,
)
from phykit.errors import PhykitUserError


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
GENE_TREES_FILE = str(SAMPLE_FILES / "gene_trees_simple.nwk")


def test_module_import_does_not_import_scipy_linalg(monkeypatch):
    module_name = "phykit.services.tree.vcv_utils"
    previous = sys.modules.pop(module_name, None)
    original_import = builtins.__import__

    def guarded_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "scipy.linalg" or name.startswith("scipy.linalg."):
            raise AssertionError(
                "vcv_utils module import should not import SciPy linalg"
            )
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", guarded_import)
    try:
        importlib.import_module(module_name)
    finally:
        imported = sys.modules.pop(module_name, None)
        if previous is not None:
            sys.modules[module_name] = previous
        parent_name, _, child_name = module_name.rpartition(".")
        parent = sys.modules.get(parent_name)
        if parent is not None:
            if previous is not None:
                setattr(parent, child_name, previous)
            elif getattr(parent, child_name, None) is imported:
                delattr(parent, child_name)


def test_module_import_does_not_import_biopython_phylo():
    code = """
import sys
import phykit.services.tree.vcv_utils as module
assert hasattr(module.Phylo, "read")
assert hasattr(module.np, "__getattr__")
assert "typing" not in sys.modules
assert "pickle" not in sys.modules
assert "numpy" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "Bio.Align" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


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

    def test_real_tree_fast_path_does_not_call_distance(self, monkeypatch):
        tree = _make_tree("(A:1.0,(B:0.5,C:0.5):0.5);")

        def fail_distance(*args, **kwargs):
            raise AssertionError("distance fallback should not be called")

        monkeypatch.setattr(tree, "distance", fail_distance)
        vcv = build_vcv_matrix(tree, ["A", "B", "C"])

        expected = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.5],
            [0.0, 0.5, 1.0],
        ])
        np.testing.assert_allclose(vcv, expected)

    def test_real_tree_fast_path_does_not_call_get_path(self, monkeypatch):
        tree = _make_tree("(A:1.0,(B:0.5,C:0.5):0.5);")

        def fail_get_path(*args, **kwargs):
            raise AssertionError("get_path should not be called")

        monkeypatch.setattr(TreeMixin, "get_path", fail_get_path)
        vcv = build_vcv_matrix(tree, ["A", "B", "C"])

        expected = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.5],
            [0.0, 0.5, 1.0],
        ])
        np.testing.assert_allclose(vcv, expected)

    def test_real_tree_fast_path_does_not_call_depths(self, monkeypatch):
        tree = _make_tree("(A:1.0,(B:0.5,C:0.5):0.5);")

        def fail_depths(*args, **kwargs):
            raise AssertionError("depths should not be called")

        monkeypatch.setattr(tree, "depths", fail_depths)
        vcv = build_vcv_matrix(tree, ["A", "B", "C"])

        expected = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.5],
            [0.0, 0.5, 1.0],
        ])
        np.testing.assert_allclose(vcv, expected)

    def test_real_tree_fast_path_does_not_call_get_terminals(self, monkeypatch):
        tree = _make_tree("(A:1.0,(B:0.5,C:0.5):0.5);")

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("get_terminals should not be called")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)
        vcv = build_vcv_matrix(tree, ["A", "B", "C"])

        expected = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.5],
            [0.0, 0.5, 1.0],
        ])
        np.testing.assert_allclose(vcv, expected)

    def test_single_tip_branches_skip_block_indexing(self, monkeypatch):
        tree = _make_tree("(A:1.0,B:2.0,C:3.0);")

        def fail_ix(*_args, **_kwargs):
            raise AssertionError("single-tip branches should update diagonals directly")

        monkeypatch.setattr(vcv_utils.np, "ix_", fail_ix)

        vcv = build_vcv_matrix(tree, ["A", "B", "C"])

        np.testing.assert_allclose(vcv, np.diag([1.0, 2.0, 3.0]))

    def test_multi_tip_branches_use_broadcast_indexing(self, monkeypatch):
        tree = _make_tree("((A:1.0,B:1.0):2.0,C:3.0);")

        def fail_ix(*_args, **_kwargs):
            raise AssertionError("multi-tip branches should not use np.ix_")

        monkeypatch.setattr(vcv_utils.np, "ix_", fail_ix)

        vcv = build_vcv_matrix(tree, ["A", "B", "C"])

        expected = np.array([
            [3.0, 2.0, 0.0],
            [2.0, 3.0, 0.0],
            [0.0, 0.0, 3.0],
        ])
        np.testing.assert_allclose(vcv, expected)

    def test_transformed_vcv_fast_path_does_not_call_get_path(self, monkeypatch):
        tree = _make_tree("(A:1.0,(B:0.5,C:0.5):0.5);")

        def fail_get_path(*args, **kwargs):
            raise AssertionError("get_path should not be called")

        monkeypatch.setattr(TreeMixin, "get_path", fail_get_path)
        vcv = build_transformed_vcv_matrix(tree, ["A", "B", "C"], lambda bl: bl * 2.0)

        expected = np.array([
            [2.0, 0.0, 0.0],
            [0.0, 2.0, 1.0],
            [0.0, 1.0, 2.0],
        ])
        np.testing.assert_allclose(vcv, expected)

    def test_transformed_vcv_fast_path_does_not_call_get_terminals(self, monkeypatch):
        tree = _make_tree("(A:1.0,(B:0.5,C:0.5):0.5);")

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("get_terminals should not be called")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)
        vcv = build_transformed_vcv_matrix(tree, ["A", "B", "C"], lambda bl: bl * 2.0)

        expected = np.array([
            [2.0, 0.0, 0.0],
            [0.0, 2.0, 1.0],
            [0.0, 1.0, 2.0],
        ])
        np.testing.assert_allclose(vcv, expected)

    def test_transformed_vcv_fast_path_matches_explicit_distance_formula(self):
        tree = _make_tree("(A:1.0,(B:0.5,C:0.5):0.5):2.0;")
        tips = ["A", "B", "C"]
        transform = lambda bl: bl * bl + 0.25

        vcv = build_transformed_vcv_matrix(tree, tips, transform)

        transformed = _make_tree("(A:1.0,(B:0.5,C:0.5):0.5):2.0;")
        for clade in transformed.find_clades():
            if clade.branch_length is not None and clade is not transformed.root:
                clade.branch_length = transform(clade.branch_length)
        expected = build_vcv_matrix(transformed, tips)

        np.testing.assert_allclose(vcv, expected)

    def test_explicit_root_branch_length_matches_distance_formula(self):
        tree = _make_tree("(A:1.0,(B:0.5,C:0.5):0.5):2.0;")
        tips = ["A", "B", "C"]
        vcv = build_vcv_matrix(tree, tips)

        expected = np.zeros((3, 3))
        root_to_tip = {name: tree.distance(tree.root, name) for name in tips}
        for i, name_i in enumerate(tips):
            expected[i, i] = root_to_tip[name_i]
            for j in range(i + 1, len(tips)):
                name_j = tips[j]
                d_ij = tree.distance(name_i, name_j)
                shared_path = (root_to_tip[name_i] + root_to_tip[name_j] - d_ij) / 2
                expected[i, j] = shared_path
                expected[j, i] = shared_path

        np.testing.assert_allclose(vcv, expected)


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

    def test_comments_blanks_and_whitespace_ignored(self, tmp_path):
        f = tmp_path / "with_comments_and_blanks.nwk"
        f.write_text(
            "   # header\n\n  (A:1.0,B:2.0,C:3.0);  \n\t# footer\n\n"
        )
        trees = parse_gene_trees(str(f))
        assert len(trees) == 1

    def test_path_list_avoids_per_row_path_objects(self, tmp_path, monkeypatch):
        (tmp_path / "one.nwk").write_text("(A:1,B:1);\n")
        (tmp_path / "two.nwk").write_text("(A:1,C:1);\n")
        gene_trees = tmp_path / "gene_tree_paths.txt"
        gene_trees.write_text("one.nwk\ntwo.nwk\n")
        real_path = Path
        parent_joins = 0
        read_paths = []

        class CountingParent:
            def __init__(self, path):
                self._path = path

            def __str__(self):
                return str(self._path)

            def __truediv__(self, other):
                nonlocal parent_joins
                parent_joins += 1
                return self._path / other

        class CountingPath:
            def __init__(self, path):
                self._path = real_path(path)

            @property
            def parent(self):
                return CountingParent(self._path.parent)

            def open(self, *args, **kwargs):
                return self._path.open(*args, **kwargs)

        def fake_read(path, fmt):
            read_paths.append(path)
            return object()

        original_import = builtins.__import__
        use_counting_path = True

        class FakePathlib:
            Path = CountingPath

        def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
            if use_counting_path and name == "pathlib" and "Path" in fromlist:
                return FakePathlib
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", fake_import)
        monkeypatch.setattr(vcv_utils.Phylo, "read", fake_read)

        try:
            trees = parse_gene_trees(str(gene_trees))
        finally:
            use_counting_path = False

        assert len(trees) == 2
        assert parent_joins == 0
        assert read_paths == [
            str(tmp_path / "one.nwk"),
            str(tmp_path / "two.nwk"),
        ]


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

    def test_shared_taxa_setup_uses_fast_tip_name_helper(self, mocker):
        from phykit.services.tree import vcv_utils

        species = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        gt = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        species_spy = mocker.spy(vcv_utils, "_get_tip_names")
        gene_spy = mocker.spy(
            vcv_utils,
            "_get_tip_names_and_missing_branch_lengths",
        )

        _, meta = build_discordance_vcv(species, [gt], ["A", "B", "C", "D"])

        assert meta["shared_taxa"] == ["A", "B", "C", "D"]
        species_spy.assert_called_once_with(species)
        gene_spy.assert_called_once_with(gt)

    def test_pruning_does_not_mutate_gene_tree(self):
        species = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        gt = _make_tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);")

        vcv, meta = build_discordance_vcv(species, [gt], ["A", "B", "C", "D"])

        assert meta["shared_taxa"] == ["A", "B", "C", "D"]
        assert {tip.name for tip in gt.get_terminals()} == {"A", "B", "C", "D", "E"}
        assert vcv.shape == (4, 4)

    def test_all_shared_gene_tree_skips_copy_prune(self, monkeypatch):
        from phykit.services.tree import vcv_utils

        species = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        gt = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")

        def fail_copy_prune(*_args, **_kwargs):
            raise AssertionError("all-shared gene trees should not be copied")

        monkeypatch.setattr(
            vcv_utils,
            "_copy_prune_gene_tree_to_shared_taxa",
            fail_copy_prune,
        )

        vcv, meta = build_discordance_vcv(species, [gt], ["A", "B", "C", "D"])

        assert meta["shared_taxa"] == ["A", "B", "C", "D"]
        assert vcv.shape == (4, 4)

    def test_pruning_uses_batch_path_for_standard_tree(self, monkeypatch):
        species = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        gt = _make_tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);")

        def fail_prune(*_args, **_kwargs):
            raise AssertionError("per-tip prune should not be used")

        monkeypatch.setattr(TreeMixin, "prune", fail_prune)

        vcv, meta = build_discordance_vcv(species, [gt], ["A", "B", "C", "D"])

        assert meta["shared_taxa"] == ["A", "B", "C", "D"]
        assert vcv.shape == (4, 4)

    def test_pruned_subset_vcv_matches_copy_prune(self):
        gt = _make_tree(
            "(((A:1,X:1):2,B:3):4,(C:1,(D:2,E:3):4):5);"
        )
        shared_ordered = ["A", "B", "C", "D"]

        pruned = _copy_prune_gene_tree_to_shared_taxa(gt, set(shared_ordered))
        expected = build_vcv_matrix(pruned, shared_ordered)
        observed = _build_pruned_subset_vcv_matrix(gt, shared_ordered)

        np.testing.assert_allclose(observed, expected)
        assert {tip.name for tip in gt.get_terminals()} == {"A", "B", "C", "D", "E", "X"}

    def test_pruned_subset_single_tip_branches_skip_block_indexing(
        self, monkeypatch
    ):
        gt = _make_tree("(A:1.0,B:2.0,X:3.0);")

        def fail_ix(*_args, **_kwargs):
            raise AssertionError("single-tip branches should update diagonals directly")

        monkeypatch.setattr(vcv_utils.np, "ix_", fail_ix)

        observed = _build_pruned_subset_vcv_matrix(gt, ["A", "B"])

        np.testing.assert_allclose(observed, np.diag([1.0, 2.0]))

    def test_pruned_subset_multi_tip_branches_use_broadcast_indexing(
        self, monkeypatch
    ):
        gt = _make_tree("((A:1.0,B:1.0):2.0,(C:1.0,X:1.0):2.0);")

        def fail_ix(*_args, **_kwargs):
            raise AssertionError("multi-tip branches should not use np.ix_")

        monkeypatch.setattr(vcv_utils.np, "ix_", fail_ix)

        observed = _build_pruned_subset_vcv_matrix(gt, ["A", "B", "C"])

        expected = np.array([
            [3.0, 2.0, 0.0],
            [2.0, 3.0, 0.0],
            [0.0, 0.0, 3.0],
        ])
        np.testing.assert_allclose(observed, expected)

    def test_extra_tip_gene_tree_skips_copy_prune_when_subset_vcv_is_available(
        self, monkeypatch
    ):
        from phykit.services.tree import vcv_utils

        species = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        gt = _make_tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);")

        def fail_copy_prune(*_args, **_kwargs):
            raise AssertionError("standard extra-tip gene trees should not be copied")

        monkeypatch.setattr(
            vcv_utils,
            "_copy_prune_gene_tree_to_shared_taxa",
            fail_copy_prune,
        )

        vcv, meta = build_discordance_vcv(species, [gt], ["A", "B", "C", "D"])

        assert meta["shared_taxa"] == ["A", "B", "C", "D"]
        assert vcv.shape == (4, 4)

    def test_branch_length_validation_uses_direct_traversal(self, monkeypatch):
        gt = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("generic traversal should not be used")

        monkeypatch.setattr(gt, "find_clades", fail_find_clades)

        assert _gene_tree_has_missing_branch_lengths(gt) is False

    def test_branch_length_validation_handles_mixed_child_counts(self, monkeypatch):
        gt = _make_tree("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);")

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("standard tree validation should use direct traversal")

        monkeypatch.setattr(gt, "find_clades", fail_find_clades)

        assert _gene_tree_has_missing_branch_lengths(gt) is False
        gt.root.clades[2].clades[1].branch_length = None
        assert _gene_tree_has_missing_branch_lengths(gt) is True

    def test_tip_names_and_branch_length_validation_share_direct_traversal(
        self, monkeypatch
    ):
        gt = _make_tree("((A:1,B):1,(C:1,D:1):1);")

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("generic traversal should not be used")

        monkeypatch.setattr(gt, "find_clades", fail_find_clades)

        tips, has_missing_branch_lengths = (
            _get_tip_names_and_missing_branch_lengths(gt)
        )

        assert tips == ["A", "B", "C", "D"]
        assert has_missing_branch_lengths is True

    def test_tip_names_and_branch_length_validation_preserves_multifurcation_order(
        self, monkeypatch
    ):
        gt = _make_tree("(A:1,(B:1,C,D:1):1,E:1);")

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("generic traversal should not be used")

        monkeypatch.setattr(gt, "find_clades", fail_find_clades)

        tips, has_missing_branch_lengths = (
            _get_tip_names_and_missing_branch_lengths(gt)
        )

        assert tips == ["A", "B", "C", "D", "E"]
        assert has_missing_branch_lengths is True

    def test_discordance_vcv_reuses_cached_branch_length_validation(
        self, monkeypatch
    ):
        species = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        gt = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")

        def fail_branch_validation(*_args, **_kwargs):
            raise AssertionError("branch validation should be cached during tip scan")

        monkeypatch.setattr(
            vcv_utils,
            "_gene_tree_has_missing_branch_lengths",
            fail_branch_validation,
        )

        vcv, meta = build_discordance_vcv(species, [gt], ["A", "B", "C", "D"])

        assert meta["shared_taxa"] == ["A", "B", "C", "D"]
        assert vcv.shape == (4, 4)

    def test_copy_prune_gene_tree_returns_original_when_all_taxa_shared(
        self, monkeypatch
    ):
        gt = _make_tree("((A:1,B:1):1,(C:1,D:1,E:1):1);")

        def fail_pickle_copy(*_args, **_kwargs):
            raise AssertionError("all-shared tree should not be copied")

        monkeypatch.setattr(vcv_utils.pickle, "dumps", fail_pickle_copy)

        pruned = _copy_prune_gene_tree_to_shared_taxa(
            gt, {"A", "B", "C", "D", "E"}
        )

        assert pruned is gt
        assert {tip.name for tip in gt.get_terminals()} == {
            "A",
            "B",
            "C",
            "D",
            "E",
        }

    def test_copy_prune_gene_tree_fallback_prunes_terminal_objects(self, monkeypatch):
        gt = _make_tree("(A:1,B:1);")
        original_prune = TreeMixin.prune

        def assert_object_prune(self, target=None, **kwargs):
            assert not isinstance(target, str)
            return original_prune(self, target=target, **kwargs)

        monkeypatch.setattr(TreeMixin, "prune", assert_object_prune)

        pruned = _copy_prune_gene_tree_to_shared_taxa(gt, {"A"})

        assert pruned is not gt
        assert {tip.name for tip in gt.get_terminals()} == {"A", "B"}
        assert {tip.name for tip in pruned.get_terminals()} == {"A"}

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

    def test_already_psd_uses_min_eigenvalue_fast_path(self, monkeypatch):
        def fail_full_eigh(*args, **kwargs):
            raise AssertionError("full eigendecomposition should not be called")

        monkeypatch.setattr(np.linalg, "eigh", fail_full_eigh)

        matrix = np.eye(4)
        result, corrected, min_eval = _nearest_psd(matrix)

        np.testing.assert_array_almost_equal(result, matrix)
        assert corrected is False
        assert min_eval == pytest.approx(1.0)

    def test_repeated_psd_checks_cache_scipy_eigvalsh_import(self, monkeypatch):
        previous_eigvalsh = vcv_utils._EIGVALSH
        vcv_utils._EIGVALSH = None
        original_import = builtins.__import__
        scipy_linalg_imports = 0

        def counting_import(name, globals=None, locals=None, fromlist=(), level=0):
            nonlocal scipy_linalg_imports
            if name == "scipy.linalg":
                scipy_linalg_imports += 1
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", counting_import)
        try:
            _nearest_psd(np.eye(4))
            _nearest_psd(np.eye(4))
        finally:
            vcv_utils._EIGVALSH = previous_eigvalsh

        assert scipy_linalg_imports == 1

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

    def test_eigendecomposition_reconstruction_matches_dense_reference_without_diag(
        self, monkeypatch
    ):
        rng = np.random.default_rng(20260623)
        Q, _ = np.linalg.qr(rng.normal(size=(12, 12)))
        eigvals = np.linspace(-0.5, 2.0, 12)
        bad = Q @ np.diag(eigvals) @ Q.T
        clipped = np.maximum(eigvals, 0.0)
        expected = Q @ np.diag(clipped) @ Q.T
        expected = (expected + expected.T) / 2.0

        original_diag = vcv_utils.np.diag

        def diag_without_vector_to_matrix(a, *args, **kwargs):
            arr = np.asarray(a)
            if arr.ndim == 1:
                raise AssertionError("PSD reconstruction should not allocate a diagonal")
            return original_diag(a, *args, **kwargs)

        monkeypatch.setattr(vcv_utils.np, "diag", diag_without_vector_to_matrix)
        monkeypatch.setattr(
            vcv_utils.np,
            "min",
            lambda *args, **kwargs: pytest.fail(
                "PSD reconstruction should use ndarray.min"
            ),
        )

        observed, corrected, min_eval = vcv_utils._nearest_psd_from_eigendecomposition(
            bad, eigvals, Q
        )

        assert corrected is True
        assert min_eval == pytest.approx(-0.5)
        np.testing.assert_allclose(observed, expected)

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
