import builtins
import importlib
import os
import json
import subprocess
import sys

import pytest
import numpy as np
from argparse import Namespace
from io import StringIO
from pathlib import Path

from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin

import phykit.services.tree.phylogenetic_ordination as phylogenetic_ordination_module
from phykit.services.tree.phylogenetic_ordination import PhylogeneticOrdination
from phykit.helpers.pgls_utils import max_lambda as compute_max_lambda
from phykit.helpers.trait_parsing import parse_multi_trait_file
from phykit.errors import PhykitUserError


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
MULTI_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_multi_traits.tsv")


def test_module_import_does_not_import_scipy_linalg_or_optimize(monkeypatch):
    module_name = "phykit.services.tree.phylogenetic_ordination"
    previous = sys.modules.pop(module_name, None)
    original_import = builtins.__import__

    def guarded_import(name, globals=None, locals=None, fromlist=(), level=0):
        if (
            name == "scipy.linalg"
            or name.startswith("scipy.linalg.")
            or name == "scipy.optimize"
            or name.startswith("scipy.optimize.")
        ):
            raise AssertionError(
                "phylogenetic_ordination module import should not import SciPy linalg/optimize"
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


def test_module_import_does_not_import_numpy_scipy_or_pgls_helper():
    code = """
import sys
import phykit.services.tree.phylogenetic_ordination as module
assert hasattr(module.np, "__getattr__")
assert callable(module.print_json)
assert callable(module.parse_multi_trait_file)
assert "pickle" not in sys.modules
assert "json" not in sys.modules
assert "typing" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
assert "phykit.helpers.trait_parsing" not in sys.modules
assert "numpy" not in sys.modules
assert "scipy.linalg" not in sys.modules
assert "scipy.optimize" not in sys.modules
assert "phykit.helpers.pgls_utils" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


@pytest.fixture
def default_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=MULTI_TRAITS_FILE,
        method="pca",
        correction="BM",
        mode="cov",
        json=False,
    )


@pytest.fixture
def corr_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=MULTI_TRAITS_FILE,
        method="pca",
        correction="BM",
        mode="corr",
        json=False,
    )


@pytest.fixture
def lambda_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=MULTI_TRAITS_FILE,
        method="pca",
        correction="lambda",
        mode="cov",
        json=False,
    )


@pytest.fixture
def lambda_corr_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=MULTI_TRAITS_FILE,
        method="pca",
        correction="lambda",
        mode="corr",
        json=False,
    )


@pytest.fixture
def tsne_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=MULTI_TRAITS_FILE,
        method="tsne",
        correction="BM",
        n_components=2,
        perplexity=None,
        n_neighbors=None,
        min_dist=0.1,
        seed=42,
        json=False,
        plot=False,
        plot_output="phylo_ordination_plot.png",
        plot_tree=False,
        color_by=None,
    )


@pytest.fixture
def umap_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=MULTI_TRAITS_FILE,
        method="umap",
        correction="BM",
        n_components=2,
        perplexity=None,
        n_neighbors=None,
        min_dist=0.1,
        seed=42,
        json=False,
        plot=False,
        plot_output="phylo_ordination_plot.png",
        plot_tree=False,
        color_by=None,
    )


class TestProcessArgs:
    def test_defaults(self):
        args = Namespace(tree="t.tre", trait_data="d.tsv")
        svc = PhylogeneticOrdination(args)
        assert svc.tree_file_path == "t.tre"
        assert svc.trait_data_path == "d.tsv"
        assert svc.method == "pca"
        assert svc.correction == "BM"
        assert svc.mode == "cov"
        assert svc.json_output is False

    def test_pca_overrides(self):
        args = Namespace(
            tree="t.tre", trait_data="d.tsv",
            json=True, method="pca", correction="lambda", mode="corr",
        )
        svc = PhylogeneticOrdination(args)
        assert svc.json_output is True
        assert svc.method == "pca"
        assert svc.correction == "lambda"
        assert svc.mode == "corr"

    def test_tsne_overrides(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            method="tsne",
            correction="lambda",
            n_components=3,
            perplexity=10.0,
            n_neighbors=5,
            min_dist=0.5,
            seed=99,
            json=True,
        )
        svc = PhylogeneticOrdination(args)
        assert svc.method == "tsne"
        assert svc.correction == "lambda"
        assert svc.n_components == 3
        assert svc.perplexity == 10.0
        assert svc.n_neighbors == 5
        assert svc.min_dist == 0.5
        assert svc.seed == 99
        assert svc.json_output is True

    def test_umap_defaults(self):
        args = Namespace(tree="t.tre", trait_data="d.tsv", method="umap")
        svc = PhylogeneticOrdination(args)
        assert svc.method == "umap"
        assert svc.correction == "BM"
        assert svc.n_components == 2
        assert svc.perplexity is None
        assert svc.n_neighbors is None
        assert svc.min_dist == 0.1
        assert svc.seed is None


class TestParseMultiTraitFile:
    def test_valid_file(self, default_args):
        svc = PhylogeneticOrdination(default_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(MULTI_TRAITS_FILE, tips)
        assert len(traits) == 8
        assert trait_names == ["body_mass", "brain_size", "longevity"]
        assert traits["raccoon"] == pytest.approx([1.04, 1.60, 2.71])
        assert traits["weasel"] == pytest.approx([-0.30, 0.85, 1.79])

    def test_missing_file(self, default_args):
        svc = PhylogeneticOrdination(default_args)
        with pytest.raises(PhykitUserError):
            parse_multi_trait_file("/nonexistent/path.tsv", ["a", "b", "c"])

    def test_non_numeric_value(self, default_args, tmp_path):
        trait_file = tmp_path / "bad.tsv"
        trait_file.write_text("taxon\ttrait1\ntaxon1\tabc\n")
        svc = PhylogeneticOrdination(default_args)
        with pytest.raises(PhykitUserError):
            parse_multi_trait_file(str(trait_file), ["taxon1"])

    def test_wrong_column_count(self, default_args, tmp_path):
        trait_file = tmp_path / "bad.tsv"
        trait_file.write_text("taxon\ttrait1\ttrait2\ntaxon1\t1.0\n")
        svc = PhylogeneticOrdination(default_args)
        with pytest.raises(PhykitUserError):
            parse_multi_trait_file(str(trait_file), ["taxon1"])

    def test_header_only(self, default_args, tmp_path):
        trait_file = tmp_path / "header_only.tsv"
        trait_file.write_text("taxon\ttrait1\ttrait2\n")
        svc = PhylogeneticOrdination(default_args)
        with pytest.raises(PhykitUserError):
            parse_multi_trait_file(str(trait_file), ["taxon1", "taxon2", "taxon3"])

    def test_comments_and_blanks(self, default_args, tmp_path):
        trait_file = tmp_path / "good.tsv"
        trait_file.write_text(
            "# comment\n\ntaxon\tt1\tt2\n"
            "taxon1\t1.0\t2.0\ntaxon2\t3.0\t4.0\ntaxon3\t5.0\t6.0\n"
        )
        svc = PhylogeneticOrdination(default_args)
        trait_names, traits = parse_multi_trait_file(
            str(trait_file), ["taxon1", "taxon2", "taxon3"]
        )
        assert len(traits) == 3
        assert trait_names == ["t1", "t2"]

    def test_taxon_mismatch_warns(self, default_args, tmp_path, capsys):
        trait_file = tmp_path / "partial.tsv"
        trait_file.write_text(
            "taxon\tt1\ntaxon1\t1.0\ntaxon2\t2.0\ntaxon3\t3.0\nextra\t4.0\n"
        )
        svc = PhylogeneticOrdination(default_args)
        trait_names, traits = parse_multi_trait_file(
            str(trait_file), ["taxon1", "taxon2", "taxon3", "taxon4"]
        )
        _, err = capsys.readouterr()
        assert "Warning" in err
        assert len(traits) == 3

    def test_too_few_shared_taxa(self, default_args, tmp_path):
        trait_file = tmp_path / "few.tsv"
        trait_file.write_text("taxon\tt1\ntaxon1\t1.0\ntaxon2\t2.0\n")
        svc = PhylogeneticOrdination(default_args)
        with pytest.raises(PhykitUserError):
            parse_multi_trait_file(str(trait_file), ["taxon1", "taxon2"])


class TestBuildVCVMatrix:
    def test_symmetric(self, default_args):
        svc = PhylogeneticOrdination(default_args)
        tree = svc.read_tree_file()
        tips = sorted(svc.get_tip_names_from_tree(tree))
        vcv = svc._build_vcv_matrix(tree, tips)
        np.testing.assert_array_almost_equal(vcv, vcv.T)

    def test_diagonal_is_root_to_tip(self, default_args):
        svc = PhylogeneticOrdination(default_args)
        tree = svc.read_tree_file()
        tips = sorted(svc.get_tip_names_from_tree(tree))
        vcv = svc._build_vcv_matrix(tree, tips)
        for i, name in enumerate(tips):
            expected = tree.distance(tree.root, name)
            assert vcv[i, i] == pytest.approx(expected, rel=1e-6)

    def test_correct_shape(self, default_args):
        svc = PhylogeneticOrdination(default_args)
        tree = svc.read_tree_file()
        tips = sorted(svc.get_tip_names_from_tree(tree))
        vcv = svc._build_vcv_matrix(tree, tips)
        assert vcv.shape == (8, 8)


class TestPhylogeneticPCACov:
    """Test BM + cov mode against R phytools::phyl.pca reference values."""

    REF_EIGENVALUES = [0.0801798810, 0.0029237493, 0.0003075388]
    REF_LOADINGS = np.array([
        [0.7347392318, -0.4226132689, 0.5306187767],
        [0.5274002700, -0.1360637256, -0.8386510703],
        [0.4266230379, 0.8960383293, 0.1229149949],
    ])
    REF_SCORES = {
        "bear":     [ 0.9781912622, -0.02928980193,  0.02678735213],
        "cat":      [-1.4420092510,  0.17646695023,  0.09307185066],
        "dog":      [-0.6517229303, -0.09142653369, -0.04614241827],
        "monkey":   [-0.6404655004,  1.09725068252, -0.20806754274],
        "raccoon":  [-0.8671208830,  0.06719921291,  0.11461079162],
        "sea_lion": [ 0.8603996846, -0.19926813823, -0.11510172228],
        "seal":     [ 0.5225838100,  0.27753903060, -0.05910744959],
        "weasel":   [-2.6397148510, -0.08880647561, -0.08051186175],
    }

    def test_eigenvalues(self, default_args):
        svc = PhylogeneticOrdination(default_args)
        tree = svc.read_tree_file()
        svc.validate_tree(tree, min_tips=3, require_branch_lengths=True)
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(MULTI_TRAITS_FILE, tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])
        vcv = svc._build_vcv_matrix(tree, ordered_names)

        C_inv = np.linalg.inv(vcv)
        ones = np.ones(n)
        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])
        Z = Y - np.outer(ones, a_hat)
        R = (Z.T @ C_inv @ Z) / (n - 1)
        eigenvalues, _ = np.linalg.eigh(R)
        eigenvalues = np.sort(eigenvalues)[::-1]

        np.testing.assert_allclose(eigenvalues, self.REF_EIGENVALUES, atol=1e-4)

    def test_loadings(self, default_args):
        svc = PhylogeneticOrdination(default_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(MULTI_TRAITS_FILE, tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])
        vcv = svc._build_vcv_matrix(tree, ordered_names)

        C_inv = np.linalg.inv(vcv)
        ones = np.ones(n)
        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])
        Z = Y - np.outer(ones, a_hat)
        R = (Z.T @ C_inv @ Z) / (n - 1)
        eigenvalues, eigenvectors = np.linalg.eigh(R)
        idx = np.argsort(eigenvalues)[::-1]
        eigenvectors = eigenvectors[:, idx]

        np.testing.assert_allclose(
            np.abs(eigenvectors), np.abs(self.REF_LOADINGS), atol=1e-4
        )

    def test_scores(self, default_args):
        svc = PhylogeneticOrdination(default_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(MULTI_TRAITS_FILE, tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])
        vcv = svc._build_vcv_matrix(tree, ordered_names)

        C_inv = np.linalg.inv(vcv)
        ones = np.ones(n)
        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])
        Z = Y - np.outer(ones, a_hat)
        R = (Z.T @ C_inv @ Z) / (n - 1)
        eigenvalues, eigenvectors = np.linalg.eigh(R)
        idx = np.argsort(eigenvalues)[::-1]
        eigenvectors = eigenvectors[:, idx]
        scores = Z @ eigenvectors

        for k, name in enumerate(ordered_names):
            ref = np.array(self.REF_SCORES[name])
            computed = scores[k]
            for pc in range(p):
                assert abs(abs(computed[pc]) - abs(ref[pc])) < 1e-4, (
                    f"Score mismatch for {name}, PC{pc+1}: "
                    f"computed={computed[pc]}, ref={ref[pc]}"
                )


class TestPhylogeneticPCACorr:
    """Test BM + corr mode against R phytools::phyl.pca reference values."""

    REF_EIGENVALUES = [2.8490059295, 0.1402470393, 0.0107470312]
    REF_LOADINGS = np.array([
        [-0.5820001250, -0.4638777874,  0.66790212815],
        [-0.5861515834, -0.3299925945, -0.73995351804],
        [-0.5636507569,  0.8221449300,  0.07984696836],
    ])
    REF_SCORES = {
        "bear":     [-5.718644630, -0.2520427517,  0.1680858429],
        "cat":      [ 8.308065874,  1.4817516853,  0.5169817433],
        "dog":      [ 4.015546433, -0.6377122605, -0.2647078216],
        "monkey":   [ 1.158323429,  7.3219170395, -1.3711746799],
        "raccoon":  [ 5.154737886,  0.7073985615,  0.6616876057],
        "sea_lion": [-4.828148632, -1.6179405137, -0.6477830531],
        "seal":     [-3.785695422,  1.7805146765, -0.3805726389],
        "weasel":   [15.786800820, -0.4888037720, -0.4806870806],
    }

    def test_eigenvalues(self, corr_args):
        svc = PhylogeneticOrdination(corr_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(MULTI_TRAITS_FILE, tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])
        vcv = svc._build_vcv_matrix(tree, ordered_names)

        C_inv = np.linalg.inv(vcv)
        ones = np.ones(n)
        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])
        Z = Y - np.outer(ones, a_hat)
        R = (Z.T @ C_inv @ Z) / (n - 1)

        D = np.sqrt(np.diag(R))
        D_inv = np.diag(1.0 / D)
        R_corr = D_inv @ R @ D_inv
        eigenvalues, _ = np.linalg.eigh(R_corr)
        eigenvalues = np.sort(eigenvalues)[::-1]

        np.testing.assert_allclose(eigenvalues, self.REF_EIGENVALUES, atol=1e-4)

    def test_loadings(self, corr_args):
        svc = PhylogeneticOrdination(corr_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(MULTI_TRAITS_FILE, tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])
        vcv = svc._build_vcv_matrix(tree, ordered_names)

        C_inv = np.linalg.inv(vcv)
        ones = np.ones(n)
        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])
        Z = Y - np.outer(ones, a_hat)
        R = (Z.T @ C_inv @ Z) / (n - 1)

        D = np.sqrt(np.diag(R))
        D_inv = np.diag(1.0 / D)
        R_corr = D_inv @ R @ D_inv
        eigenvalues, eigenvectors = np.linalg.eigh(R_corr)
        idx = np.argsort(eigenvalues)[::-1]
        eigenvectors = eigenvectors[:, idx]

        np.testing.assert_allclose(
            np.abs(eigenvectors), np.abs(self.REF_LOADINGS), atol=1e-4
        )

    def test_scores(self, corr_args):
        svc = PhylogeneticOrdination(corr_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(MULTI_TRAITS_FILE, tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])
        vcv = svc._build_vcv_matrix(tree, ordered_names)

        C_inv = np.linalg.inv(vcv)
        ones = np.ones(n)
        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])
        Z = Y - np.outer(ones, a_hat)
        R = (Z.T @ C_inv @ Z) / (n - 1)

        D = np.sqrt(np.diag(R))
        D_inv = np.diag(1.0 / D)
        R_corr = D_inv @ R @ D_inv
        eigenvalues, eigenvectors = np.linalg.eigh(R_corr)
        idx = np.argsort(eigenvalues)[::-1]
        eigenvectors = eigenvectors[:, idx]

        Z_std = Z @ D_inv
        scores = Z_std @ eigenvectors

        for k, name in enumerate(ordered_names):
            ref = np.array(self.REF_SCORES[name])
            computed = scores[k]
            for pc in range(p):
                assert abs(abs(computed[pc]) - abs(ref[pc])) < 1e-4, (
                    f"Score mismatch for {name}, PC{pc+1}: "
                    f"computed={computed[pc]}, ref={ref[pc]}"
                )

    def test_run_pca_corr_mode_avoids_dense_diagonal_scaling(
        self, corr_args, monkeypatch, capsys
    ):
        import phykit.services.tree.phylogenetic_ordination as po

        svc = PhylogeneticOrdination(corr_args)
        svc.json_output = True
        rng = np.random.default_rng(20260623)
        n = 12
        p = 5
        Z = rng.normal(size=(n, p))
        C_inv_Z = Z.copy()

        original_diag = po.np.diag

        def diag_without_vector_to_matrix(a, *args, **kwargs):
            arr = np.asarray(a)
            if arr.ndim == 1:
                raise AssertionError("corr-mode PCA should not allocate D_inv")
            return original_diag(a, *args, **kwargs)

        monkeypatch.setattr(po.np, "diag", diag_without_vector_to_matrix)

        svc._run_pca(
            Z,
            C_inv_Z,
            n,
            p,
            tree=None,
            ordered_names=[f"taxon{i}" for i in range(n)],
            trait_names=[f"trait{i}" for i in range(p)],
            Y=Z,
            lambda_val=None,
            log_likelihood=None,
        )

        payload = json.loads(capsys.readouterr().out)
        assert set(payload["scores"]) == {f"taxon{i}" for i in range(n)}


class TestPhylogeneticPCALambda:
    """Test lambda + cov mode against R phytools::phyl.pca reference values."""

    REF_LAMBDA = 6.610696e-05
    REF_LOG_LIKELIHOOD = -6.016939
    REF_EIGENVALUES = [0.076301214816, 0.002241847818, 0.000315397635]

    def test_multi_trait_log_likelihood_cholesky_matches_inverse(self, lambda_args):
        rng = np.random.default_rng(20260623)
        svc = PhylogeneticOrdination(lambda_args)
        n = 16
        p = 4
        A = rng.normal(size=(n, n))
        C = A @ A.T + np.eye(n)
        Y = rng.normal(size=(n, p))

        fast = svc._multi_trait_log_likelihood_cholesky(Y, C)
        inverse = svc._multi_trait_log_likelihood_inverse(Y, C)

        assert fast == pytest.approx(inverse)

    def test_multi_trait_log_likelihood_cholesky_uses_single_solve(
        self, lambda_args, monkeypatch
    ):
        rng = np.random.default_rng(20260628)
        svc = PhylogeneticOrdination(lambda_args)
        n = 16
        p = 4
        A = rng.normal(size=(n, n))
        C = A @ A.T + np.eye(n)
        Y = rng.normal(size=(n, p))
        original = phylogenetic_ordination_module.cho_solve
        calls = 0

        def counting_cho_solve(*args, **kwargs):
            nonlocal calls
            calls += 1
            return original(*args, **kwargs)

        monkeypatch.setattr(
            phylogenetic_ordination_module, "cho_solve", counting_cho_solve
        )

        fast = svc._multi_trait_log_likelihood_cholesky(Y, C)
        inverse = svc._multi_trait_log_likelihood_inverse(Y, C)

        assert calls == 1
        assert fast == pytest.approx(inverse)

    def test_center_traits_cholesky_matches_inverse(self, lambda_args):
        rng = np.random.default_rng(20260623)
        svc = PhylogeneticOrdination(lambda_args)
        n = 16
        p = 4
        A = rng.normal(size=(n, n))
        C = A @ A.T + np.eye(n)
        Y = rng.normal(size=(n, p))

        Z_fast, C_inv_Z_fast = svc._center_traits_by_vcv_cholesky(
            Y, C, include_weighted=True
        )
        Z_inv, C_inv_Z_inv = svc._center_traits_by_vcv_inverse(
            Y, C, include_weighted=True
        )

        np.testing.assert_allclose(Z_fast, Z_inv)
        np.testing.assert_allclose(C_inv_Z_fast, C_inv_Z_inv)

    def test_center_traits_cholesky_uses_single_solve(self, lambda_args, monkeypatch):
        rng = np.random.default_rng(20260628)
        svc = PhylogeneticOrdination(lambda_args)
        n = 16
        p = 4
        A = rng.normal(size=(n, n))
        C = A @ A.T + np.eye(n)
        Y = rng.normal(size=(n, p))
        original = phylogenetic_ordination_module.cho_solve
        calls = 0

        def counting_cho_solve(*args, **kwargs):
            nonlocal calls
            calls += 1
            return original(*args, **kwargs)

        monkeypatch.setattr(
            phylogenetic_ordination_module, "cho_solve", counting_cho_solve
        )

        Z_fast, C_inv_Z_fast = svc._center_traits_by_vcv_cholesky(
            Y, C, include_weighted=True
        )
        Z_inv, C_inv_Z_inv = svc._center_traits_by_vcv_inverse(
            Y, C, include_weighted=True
        )

        assert calls == 1
        np.testing.assert_allclose(Z_fast, Z_inv)
        np.testing.assert_allclose(C_inv_Z_fast, C_inv_Z_inv)

    def test_center_traits_inverse_matches_scalar_reference(self, lambda_args):
        import phykit.services.tree.phylogenetic_ordination as po

        rng = np.random.default_rng(20260623)
        svc = PhylogeneticOrdination(lambda_args)
        n = 9
        p = 7
        A = rng.normal(size=(n, n))
        C = A @ A.T + np.eye(n)
        Y = rng.normal(size=(n, p))

        C_inv = np.linalg.inv(C)
        ones = np.ones(n)
        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])
        expected_Z = Y - np.outer(ones, a_hat)
        expected_C_inv_Z = C_inv @ expected_Z

        def fail_outer(*args, **kwargs):
            raise AssertionError("centering should use broadcasting, not np.outer")

        original_outer = po.np.outer
        po.np.outer = fail_outer
        try:
            observed_Z, observed_C_inv_Z = svc._center_traits_by_vcv_inverse(
                Y, C, include_weighted=True
            )
        finally:
            po.np.outer = original_outer

        np.testing.assert_allclose(observed_Z, expected_Z)
        np.testing.assert_allclose(observed_C_inv_Z, expected_C_inv_Z)

    def test_center_traits_cholesky_uses_broadcast_centering(
        self, lambda_args, monkeypatch
    ):
        import phykit.services.tree.phylogenetic_ordination as po

        rng = np.random.default_rng(20260623)
        svc = PhylogeneticOrdination(lambda_args)
        n = 12
        p = 6
        A = rng.normal(size=(n, n))
        C = A @ A.T + np.eye(n)
        Y = rng.normal(size=(n, p))

        def fail_outer(*args, **kwargs):
            raise AssertionError("centering should use broadcasting, not np.outer")

        monkeypatch.setattr(po.np, "outer", fail_outer)
        observed_Z, observed_C_inv_Z = svc._center_traits_by_vcv_inverse(
            Y, C, include_weighted=True
        )

        assert observed_Z.shape == Y.shape
        assert observed_C_inv_Z.shape == Y.shape
        svc._center_traits_by_vcv_cholesky(Y, C, include_weighted=True)

    def test_multi_trait_log_likelihood_inverse_matches_scalar_reference(
        self, lambda_args
    ):
        rng = np.random.default_rng(20260623)
        svc = PhylogeneticOrdination(lambda_args)
        n = 9
        p = 7
        A = rng.normal(size=(n, n))
        C = A @ A.T + np.eye(n)
        Y = rng.normal(size=(n, p))

        C_inv = np.linalg.inv(C)
        ones = np.ones(n)
        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])
        Z = Y - np.outer(ones, a_hat)
        R_mle = (Z.T @ C_inv @ Z) / n
        sign_C, logdet_C = np.linalg.slogdet(C)
        sign_R, logdet_R = np.linalg.slogdet(R_mle)
        expected = -0.5 * (
            n * p * np.log(2 * np.pi) + p * logdet_C + n * logdet_R + n * p
        )

        assert sign_C > 0
        assert sign_R > 0
        assert svc._multi_trait_log_likelihood_inverse(Y, C) == pytest.approx(
            expected
        )

    def test_lambda_value(self, lambda_args):
        svc = PhylogeneticOrdination(lambda_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(MULTI_TRAITS_FILE, tips)

        ordered_names = sorted(traits.keys())
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])
        vcv = svc._build_vcv_matrix(tree, ordered_names)
        max_lam = compute_max_lambda(tree)

        lambda_val, log_likelihood = svc._multi_trait_lambda(Y, vcv, max_lam)

        assert lambda_val == pytest.approx(self.REF_LAMBDA, abs=1e-3)
        assert log_likelihood == pytest.approx(self.REF_LOG_LIKELIHOOD, abs=0.2)

    def test_eigenvalues(self, lambda_args):
        svc = PhylogeneticOrdination(lambda_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(MULTI_TRAITS_FILE, tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])
        vcv = svc._build_vcv_matrix(tree, ordered_names)
        max_lam = compute_max_lambda(tree)

        lambda_val, _ = svc._multi_trait_lambda(Y, vcv, max_lam)

        diag_vals = np.diag(vcv).copy()
        vcv_t = vcv * lambda_val
        np.fill_diagonal(vcv_t, diag_vals)

        C_inv = np.linalg.inv(vcv_t)
        ones = np.ones(n)
        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])
        Z = Y - np.outer(ones, a_hat)
        R = (Z.T @ C_inv @ Z) / (n - 1)

        eigenvalues, _ = np.linalg.eigh(R)
        eigenvalues = np.sort(eigenvalues)[::-1]

        np.testing.assert_allclose(eigenvalues, self.REF_EIGENVALUES, atol=1e-4)


class TestTSNEEmbedding:
    def test_shape(self, tsne_args):
        svc = PhylogeneticOrdination(tsne_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(MULTI_TRAITS_FILE, tips)
        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])
        vcv = svc._build_vcv_matrix(tree, ordered_names)
        C_inv = np.linalg.inv(vcv)
        ones = np.ones(n)
        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])
        Z = Y - np.outer(ones, a_hat)

        embedding, params = svc._embed_tsne(Z, n)
        assert embedding.shape == (n, 2)

    def test_deterministic_with_seed(self, tsne_args):
        svc = PhylogeneticOrdination(tsne_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(MULTI_TRAITS_FILE, tips)
        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])
        vcv = svc._build_vcv_matrix(tree, ordered_names)
        C_inv = np.linalg.inv(vcv)
        ones = np.ones(n)
        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])
        Z = Y - np.outer(ones, a_hat)

        embedding1, _ = svc._embed_tsne(Z, n)
        embedding2, _ = svc._embed_tsne(Z, n)
        np.testing.assert_array_almost_equal(embedding1, embedding2)

    def test_auto_perplexity(self, tsne_args):
        svc = PhylogeneticOrdination(tsne_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(MULTI_TRAITS_FILE, tips)
        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])
        vcv = svc._build_vcv_matrix(tree, ordered_names)
        C_inv = np.linalg.inv(vcv)
        ones = np.ones(n)
        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])
        Z = Y - np.outer(ones, a_hat)

        _, params = svc._embed_tsne(Z, n)
        expected_perplexity = min(30.0, (n - 1) / 3.0)
        assert params["perplexity"] == round(expected_perplexity, 2)


class TestUMAPEmbedding:
    def test_shape(self, umap_args):
        svc = PhylogeneticOrdination(umap_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(MULTI_TRAITS_FILE, tips)
        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])
        vcv = svc._build_vcv_matrix(tree, ordered_names)
        C_inv = np.linalg.inv(vcv)
        ones = np.ones(n)
        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])
        Z = Y - np.outer(ones, a_hat)

        embedding, params = svc._embed_umap(Z, n)
        assert embedding.shape == (n, 2)

    def test_deterministic_with_seed(self, umap_args):
        svc = PhylogeneticOrdination(umap_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(MULTI_TRAITS_FILE, tips)
        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])
        vcv = svc._build_vcv_matrix(tree, ordered_names)
        C_inv = np.linalg.inv(vcv)
        ones = np.ones(n)
        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])
        Z = Y - np.outer(ones, a_hat)

        embedding1, _ = svc._embed_umap(Z, n)
        embedding2, _ = svc._embed_umap(Z, n)
        np.testing.assert_array_almost_equal(embedding1, embedding2)

    def test_auto_n_neighbors(self, umap_args):
        svc = PhylogeneticOrdination(umap_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(MULTI_TRAITS_FILE, tips)
        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])
        vcv = svc._build_vcv_matrix(tree, ordered_names)
        C_inv = np.linalg.inv(vcv)
        ones = np.ones(n)
        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])
        Z = Y - np.outer(ones, a_hat)

        _, params = svc._embed_umap(Z, n)
        expected_neighbors = min(15, n - 1)
        assert params["n_neighbors"] == expected_neighbors


class TestLambdaCorrection:
    def test_lambda_produces_different_embedding(self, tsne_args):
        svc_bm = PhylogeneticOrdination(tsne_args)
        tree = svc_bm.read_tree_file()
        tips = svc_bm.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(MULTI_TRAITS_FILE, tips)
        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])

        # BM centering
        vcv_bm = svc_bm._build_vcv_matrix(tree, ordered_names)
        C_inv_bm = np.linalg.inv(vcv_bm)
        ones = np.ones(n)
        denom_bm = ones @ C_inv_bm @ ones
        a_hat_bm = np.array([(ones @ C_inv_bm @ Y[:, j]) / denom_bm for j in range(p)])
        Z_bm = Y - np.outer(ones, a_hat_bm)

        # Lambda centering
        vcv_lam = svc_bm._build_vcv_matrix(tree, ordered_names)
        max_lam = compute_max_lambda(tree)
        lambda_val, _ = svc_bm._multi_trait_lambda(Y, vcv_lam, max_lam)
        diag_vals = np.diag(vcv_lam).copy()
        vcv_lam = vcv_lam * lambda_val
        np.fill_diagonal(vcv_lam, diag_vals)
        C_inv_lam = np.linalg.inv(vcv_lam)
        denom_lam = ones @ C_inv_lam @ ones
        a_hat_lam = np.array([(ones @ C_inv_lam @ Y[:, j]) / denom_lam for j in range(p)])
        Z_lam = Y - np.outer(ones, a_hat_lam)

        assert not np.allclose(Z_bm, Z_lam)

    def test_lambda_value_returned(self):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="tsne",
            correction="lambda",
            n_components=2,
            perplexity=None,
            n_neighbors=None,
            min_dist=0.1,
            seed=42,
            json=True,
        )
        svc = PhylogeneticOrdination(args)
        svc.run()


class TestSmallSampleGuards:
    def test_tsne_too_few_taxa(self, tmp_path):
        tree_file = tmp_path / "small.tre"
        tree_file.write_text("((A:1,B:1):1,C:2);")
        trait_file = tmp_path / "small.tsv"
        trait_file.write_text("taxon\tt1\tt2\nA\t1.0\t2.0\nB\t3.0\t4.0\nC\t5.0\t6.0\n")
        args = Namespace(
            tree=str(tree_file),
            trait_data=str(trait_file),
            method="tsne",
            correction="BM",
            n_components=2,
            perplexity=None,
            n_neighbors=None,
            min_dist=0.1,
            seed=42,
            json=False,
            plot=False,
            plot_output="test.png",
            plot_tree=False,
            color_by=None,
        )
        svc = PhylogeneticOrdination(args)
        with pytest.raises(PhykitUserError) as exc_info:
            svc.run()
        assert any("perplexity" in m for m in exc_info.value.messages)

    def test_umap_too_few_taxa(self, tmp_path):
        tree_file = tmp_path / "tiny.tre"
        tree_file.write_text("(A:1,B:1);")
        trait_file = tmp_path / "tiny.tsv"
        trait_file.write_text("taxon\tt1\tt2\nA\t1.0\t2.0\nB\t3.0\t4.0\n")
        args = Namespace(
            tree=str(tree_file),
            trait_data=str(trait_file),
            method="umap",
            correction="BM",
            n_components=2,
            perplexity=None,
            n_neighbors=None,
            min_dist=0.1,
            seed=42,
            json=False,
            plot=False,
            plot_output="test.png",
            plot_tree=False,
            color_by=None,
        )
        svc = PhylogeneticOrdination(args)
        with pytest.raises(PhykitUserError):
            svc.run()


class TestTreeColorTrait:
    def test_resolve_tree_color_trait_streams_color_file_once(self, monkeypatch):
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
                        "A\t1.0\textra\n",
                        "B\t2.0\n",
                        "C\t3.0\n",
                    ]
                )

        open_calls = 0

        def fake_open(path):
            nonlocal open_calls
            assert path == "colors.tsv"
            open_calls += 1
            return StreamingOnlyFile()

        svc = PhylogeneticOrdination.__new__(PhylogeneticOrdination)
        captured = {}

        def reconstruct(tree, trait_col, taxon_names):
            captured["trait_col"] = trait_col.copy()
            captured["taxon_names"] = list(taxon_names)
            return {0: np.array([trait_col[0, 0]])}, {}, tree

        monkeypatch.setattr(os.path, "isfile", lambda path: path == "colors.tsv")
        monkeypatch.setattr("builtins.open", fake_open)
        monkeypatch.setattr(svc, "_reconstruct_ancestral_scores", reconstruct)

        anc_vals, node_distances, tree_pruned = svc._resolve_tree_color_trait(
            "colors.tsv",
            ["x", "y"],
            np.zeros((3, 2)),
            ["B", "A", "C"],
            "tree",
        )

        assert open_calls == 1
        assert captured["taxon_names"] == ["B", "A", "C"]
        np.testing.assert_allclose(captured["trait_col"].ravel(), [2.0, 1.0, 3.0])
        assert anc_vals == {0: 2.0}
        assert node_distances == {}
        assert tree_pruned == "tree"


class TestRun:
    def test_run_uses_unmodified_tree_read(self, mocker):
        args = Namespace(
            tree="/some/path/to/file.tre",
            trait_data="/some/path/to/traits.tsv",
            method="pca",
            correction="BM",
            mode="cov",
            n_components=2,
            perplexity=None,
            n_neighbors=None,
            min_dist=0.1,
            seed=None,
            json=False,
            plot=False,
            plot_output="plot.png",
            plot_tree=False,
            no_plot_tree=True,
            color_by=None,
            tree_color_by=None,
            gene_trees=None,
        )
        svc = PhylogeneticOrdination(args)
        tree = object()
        mocked_read = mocker.patch.object(
            PhylogeneticOrdination,
            "read_tree_file_unmodified",
            return_value=tree,
        )
        mocker.patch.object(
            PhylogeneticOrdination,
            "read_tree_file",
            side_effect=AssertionError("run should not copy cached trees"),
        )
        mocked_validate = mocker.patch.object(
            PhylogeneticOrdination, "validate_tree"
        )
        mocker.patch.object(
            PhylogeneticOrdination,
            "get_tip_names_from_tree",
            return_value=["a", "b", "c"],
        )
        mocker.patch(
            "phykit.services.tree.phylogenetic_ordination.parse_multi_trait_file",
            return_value=(
                ["x", "y"],
                {
                    "a": [1.0, 2.0],
                    "b": [2.0, 3.0],
                    "c": [3.0, 4.0],
                },
            ),
        )
        mocked_vcv = mocker.patch(
            "phykit.services.tree.vcv_utils.build_vcv_matrix",
            return_value=np.eye(3),
        )
        mocker.patch.object(
            PhylogeneticOrdination,
            "_center_traits_by_vcv",
            return_value=(np.zeros((3, 2)), np.zeros((3, 2))),
        )
        mocked_pca = mocker.patch.object(PhylogeneticOrdination, "_run_pca")

        svc.run()

        mocked_read.assert_called_once_with()
        mocked_validate.assert_called_once_with(
            tree,
            min_tips=3,
            require_branch_lengths=True,
            context="phylogenetic ordination",
        )
        mocked_vcv.assert_called_once_with(tree, ["a", "b", "c"])
        assert mocked_pca.call_args.args[4] is tree

    def test_pca_cov_text_output(self, default_args, capsys):
        svc = PhylogeneticOrdination(default_args)
        svc.run()
        out, _ = capsys.readouterr()
        assert "Eigenvalues:" in out
        assert "Loadings:" in out
        assert "Scores:" in out
        assert "PC1" in out
        assert "body_mass" in out
        assert "raccoon" in out

    def test_pca_corr_text_output(self, corr_args, capsys):
        svc = PhylogeneticOrdination(corr_args)
        svc.run()
        out, _ = capsys.readouterr()
        assert "Eigenvalues:" in out
        assert "Loadings:" in out
        assert "Scores:" in out

    def test_pca_lambda_text_output(self, lambda_args, capsys):
        svc = PhylogeneticOrdination(lambda_args)
        svc.run()
        out, _ = capsys.readouterr()
        assert "Lambda:" in out
        assert "Log-likelihood:" in out

    def test_pca_text_output_batches_rows(self, default_args, mocker):
        svc = PhylogeneticOrdination(default_args)
        printed = mocker.patch("builtins.print")

        svc._print_pca_text_output(
            eigenvalues=np.array([2.0, 1.0]),
            proportions=np.array([0.6666667, 0.3333333]),
            eigenvectors=np.array([[0.1, 0.2], [0.3, 0.4]]),
            scores=np.array([[1.2, 1.3], [2.2, 2.3]]),
            trait_names=["trait_a", "trait_b"],
            taxon_names=["taxon_a", "taxon_b"],
            pc_labels=["PC1", "PC2"],
            lambda_val=0.5,
            log_likelihood=-12.25,
        )

        printed.assert_called_once_with(
            "Eigenvalues:\n"
            "\tPC1\tPC2\n"
            "eigenvalue\t2.000000\t1.000000\n"
            "proportion\t0.666667\t0.333333\n"
            "\n"
            "Loadings:\n"
            "\tPC1\tPC2\n"
            "trait_a\t0.100000\t0.200000\n"
            "trait_b\t0.300000\t0.400000\n"
            "\n"
            "Scores:\n"
            "\tPC1\tPC2\n"
            "taxon_a\t1.200000\t1.300000\n"
            "taxon_b\t2.200000\t2.300000\n"
            "\n"
            "Lambda: 0.500000\n"
            "Log-likelihood: -12.250000"
        )

    def test_format_pca_result_preserves_payload_shape(self, default_args):
        svc = PhylogeneticOrdination(default_args)

        result = svc._format_pca_result(
            eigenvalues=np.array([2.0, 1.0]),
            proportions=np.array([0.6666667, 0.3333333]),
            eigenvectors=np.array([[0.1, 0.2], [0.3, 0.4]]),
            scores=np.array([[1.2, 1.3], [2.2, 2.3]]),
            trait_names=["trait_a", "trait_b"],
            taxon_names=["taxon_a", "taxon_b"],
            pc_labels=["PC1", "PC2"],
            lambda_val=0.5,
            log_likelihood=-12.25,
        )

        assert result == {
            "eigenvalues": {"PC1": 2.0, "PC2": 1.0},
            "proportion_of_variance": {
                "PC1": 0.6666667,
                "PC2": 0.3333333,
            },
            "loadings": {
                "trait_a": {"PC1": 0.1, "PC2": 0.2},
                "trait_b": {"PC1": 0.3, "PC2": 0.4},
            },
            "scores": {
                "taxon_a": {"PC1": 1.2, "PC2": 1.3},
                "taxon_b": {"PC1": 2.2, "PC2": 2.3},
            },
            "lambda": 0.5,
            "log_likelihood": -12.25,
        }

    def test_pca_json_output(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="pca",
            correction="BM",
            mode="cov",
            json=True,
        )
        svc = PhylogeneticOrdination(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "eigenvalues" in data
        assert "proportion_of_variance" in data
        assert "loadings" in data
        assert "scores" in data
        assert "lambda" not in data

    def test_pca_lambda_json_output(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="pca",
            correction="lambda",
            mode="cov",
            json=True,
        )
        svc = PhylogeneticOrdination(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "eigenvalues" in data
        assert "lambda" in data
        assert "log_likelihood" in data

    def test_pca_eigenvalues_match_reference(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="pca",
            correction="BM",
            mode="cov",
            json=True,
        )
        svc = PhylogeneticOrdination(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        ev = data["eigenvalues"]
        assert ev["PC1"] == pytest.approx(0.0801798810, abs=1e-4)
        assert ev["PC2"] == pytest.approx(0.0029237493, abs=1e-4)
        assert ev["PC3"] == pytest.approx(0.0003075388, abs=1e-4)

    def test_tsne_text_output(self, tsne_args, capsys):
        svc = PhylogeneticOrdination(tsne_args)
        svc.run()
        out, _ = capsys.readouterr()
        assert "Method: tsne" in out
        assert "Correction: BM" in out
        assert "Perplexity:" in out
        assert "Embedding:" in out
        assert "Dim1" in out
        assert "Dim2" in out
        assert "raccoon" in out

    def test_umap_text_output(self, umap_args, capsys):
        svc = PhylogeneticOrdination(umap_args)
        svc.run()
        out, _ = capsys.readouterr()
        assert "Method: umap" in out
        assert "n_neighbors:" in out
        assert "Embedding:" in out

    def test_dimreduce_text_output_batches_rows(self, tsne_args, mocker):
        svc = PhylogeneticOrdination(tsne_args)
        printed = mocker.patch("builtins.print")

        svc._print_dimreduce_text_output(
            embedding=np.array([[1.2, 1.3], [2.2, 2.3]]),
            taxon_names=["taxon_a", "taxon_b"],
            params={"perplexity": 3.5, "min_dist": 0.2},
            lambda_val=0.5,
            log_likelihood=-12.25,
        )

        printed.assert_called_once_with(
            "Method: tsne\n"
            "Correction: BM\n"
            "Perplexity: 3.5\n"
            "min_dist: 0.2\n"
            "\n"
            "Embedding:\n"
            "\tDim1\tDim2\n"
            "taxon_a\t1.200000\t1.300000\n"
            "taxon_b\t2.200000\t2.300000\n"
            "\n"
            "Lambda: 0.500000\n"
            "Log-likelihood: -12.250000"
        )

    def test_tsne_json_output(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="tsne",
            correction="BM",
            n_components=2,
            perplexity=None,
            n_neighbors=None,
            min_dist=0.1,
            seed=42,
            json=True,
        )
        svc = PhylogeneticOrdination(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert data["method"] == "tsne"
        assert data["correction"] == "BM"
        assert "embedding" in data
        assert "bear" in data["embedding"]
        assert "Dim1" in data["embedding"]["bear"]
        assert "Dim2" in data["embedding"]["bear"]
        assert "lambda" not in data

    def test_dimreduce_lambda_json_output(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="tsne",
            correction="lambda",
            n_components=2,
            perplexity=None,
            n_neighbors=None,
            min_dist=0.1,
            seed=42,
            json=True,
        )
        svc = PhylogeneticOrdination(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "lambda" in data
        assert "log_likelihood" in data
        assert "embedding" in data

    def test_subset_check(self, tsne_args, capsys):
        svc = PhylogeneticOrdination(tsne_args)
        svc.run()
        out, _ = capsys.readouterr()
        for taxon in ["bear", "cat", "dog", "monkey", "raccoon", "sea_lion", "seal", "weasel"]:
            assert taxon in out

    def test_no_plot_by_default(self, tmp_path, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="pca",
            correction="BM",
            mode="cov",
            json=False,
        )
        svc = PhylogeneticOrdination(args)
        svc.run()
        out, _ = capsys.readouterr()
        assert "Saved PCA plot:" not in out


class TestPlot:
    def test_pca_plot_creates_file(self, tmp_path, capsys):
        plot_path = str(tmp_path / "test_pca.png")
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="pca",
            correction="BM",
            mode="cov",
            json=False,
            plot=True,
            plot_output=plot_path,
            plot_tree=False,
            color_by=None,
        )
        svc = PhylogeneticOrdination(args)
        svc.run()
        assert os.path.exists(plot_path)
        assert os.path.getsize(plot_path) > 0
        out, _ = capsys.readouterr()
        assert "Saved PCA plot:" in out

    def test_pca_plot_json_includes_plot_output(self, tmp_path, capsys):
        plot_path = str(tmp_path / "test_pca.png")
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="pca",
            correction="BM",
            mode="cov",
            json=True,
            plot=True,
            plot_output=plot_path,
            plot_tree=False,
            color_by=None,
        )
        svc = PhylogeneticOrdination(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "plot_output" in data
        assert data["plot_output"] == plot_path

    def test_pca_plot_tree_creates_file(self, tmp_path, capsys):
        plot_path = str(tmp_path / "test_tree.png")
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="pca",
            correction="BM",
            mode="cov",
            json=False,
            plot=True,
            plot_output=plot_path,
            plot_tree=True,
            color_by=None,
        )
        svc = PhylogeneticOrdination(args)
        svc.run()
        assert os.path.exists(plot_path)
        assert os.path.getsize(plot_path) > 0

    def test_pca_plot_color_by_column(self, tmp_path, capsys):
        plot_path = str(tmp_path / "test_color.png")
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="pca",
            correction="BM",
            mode="cov",
            json=False,
            plot=True,
            plot_output=plot_path,
            plot_tree=False,
            color_by="body_mass",
        )
        svc = PhylogeneticOrdination(args)
        svc.run()
        assert os.path.exists(plot_path)
        assert os.path.getsize(plot_path) > 0

    def test_pca_plot_tree_and_color_by_together(self, tmp_path, capsys):
        plot_path = str(tmp_path / "test_both.png")
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="pca",
            correction="BM",
            mode="cov",
            json=False,
            plot=True,
            plot_output=plot_path,
            plot_tree=True,
            color_by="body_mass",
        )
        svc = PhylogeneticOrdination(args)
        svc.run()
        assert os.path.exists(plot_path)
        assert os.path.getsize(plot_path) > 0

    def test_dimreduce_plot_created(self, tmp_path, capsys):
        plot_path = str(tmp_path / "test_dimreduce.png")
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="tsne",
            correction="BM",
            n_components=2,
            perplexity=None,
            n_neighbors=None,
            min_dist=0.1,
            seed=42,
            json=False,
            plot=True,
            plot_output=plot_path,
            plot_tree=False,
            color_by=None,
        )
        svc = PhylogeneticOrdination(args)
        svc.run()
        assert os.path.exists(plot_path)
        assert os.path.getsize(plot_path) > 0
        out, _ = capsys.readouterr()
        assert "Saved plot:" in out

    def test_dimreduce_plot_tree_created(self, tmp_path, capsys):
        plot_path = str(tmp_path / "test_tree.png")
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="tsne",
            correction="BM",
            n_components=2,
            perplexity=None,
            n_neighbors=None,
            min_dist=0.1,
            seed=42,
            json=False,
            plot=True,
            plot_output=plot_path,
            plot_tree=True,
            color_by=None,
        )
        svc = PhylogeneticOrdination(args)
        svc.run()
        assert os.path.exists(plot_path)
        assert os.path.getsize(plot_path) > 0

    def test_dimreduce_color_by_works(self, tmp_path, capsys):
        plot_path = str(tmp_path / "test_color.png")
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="tsne",
            correction="BM",
            n_components=2,
            perplexity=None,
            n_neighbors=None,
            min_dist=0.1,
            seed=42,
            json=False,
            plot=True,
            plot_output=plot_path,
            plot_tree=False,
            color_by="body_mass",
        )
        svc = PhylogeneticOrdination(args)
        svc.run()
        assert os.path.exists(plot_path)
        assert os.path.getsize(plot_path) > 0


class TestReconstructAncestralScores:
    def test_prune_setup_uses_fast_tip_names_and_set_membership(
        self, default_args, mocker
    ):
        class ContainsFailList(list):
            def __contains__(self, item):
                raise AssertionError("ordered_names should be converted to a set")

        svc = PhylogeneticOrdination(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        tree_copy = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        ordered_names = ContainsFailList(["A", "C", "D"])
        scores = np.array([[0.0, 0.0], [2.0, 2.0], [4.0, 4.0]])
        spy = mocker.spy(svc, "get_tip_names_from_tree")
        fast_copy = mocker.patch.object(svc, "_fast_copy", return_value=tree_copy)
        prune = mocker.spy(svc, "prune_tree_using_taxa_list")

        _, _, tree_pruned = svc._reconstruct_ancestral_scores(
            tree, scores, ordered_names
        )

        assert spy.call_count == 1
        fast_copy.assert_called_once_with(tree)
        prune.assert_called_once_with(tree_copy, ["B"])
        assert {tip.name for tip in tree_pruned.get_terminals()} == {"A", "C", "D"}

    def test_reconstruct_skips_copy_when_all_tips_are_shared(
        self, default_args, mocker
    ):
        svc = PhylogeneticOrdination(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        ordered_names = ["A", "B", "C", "D"]
        scores = np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 2.0], [3.0, 3.0]])
        fast_copy = mocker.patch.object(
            svc,
            "_fast_copy",
            side_effect=AssertionError("all-shared reconstruction should not copy tree"),
        )

        _, _, tree_pruned = svc._reconstruct_ancestral_scores(
            tree, scores, ordered_names
        )

        fast_copy.assert_not_called()
        assert tree_pruned is tree

    def test_all_nodes_get_scores(self, default_args):
        svc = PhylogeneticOrdination(default_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(MULTI_TRAITS_FILE, tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])
        vcv = svc._build_vcv_matrix(tree, ordered_names)

        C_inv = np.linalg.inv(vcv)
        ones = np.ones(n)
        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])
        Z = Y - np.outer(ones, a_hat)
        R = (Z.T @ C_inv @ Z) / (n - 1)
        eigenvalues, eigenvectors = np.linalg.eigh(R)
        idx = np.argsort(eigenvalues)[::-1]
        eigenvectors = eigenvectors[:, idx]
        scores = Z @ eigenvectors

        node_estimates, node_distances, tree_pruned = \
            svc._reconstruct_ancestral_scores(tree, scores, ordered_names)

        all_clades = list(tree_pruned.find_clades())
        for clade in all_clades:
            assert id(clade) in node_estimates

    def test_tips_match_original_scores(self, default_args):
        svc = PhylogeneticOrdination(default_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(MULTI_TRAITS_FILE, tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])
        vcv = svc._build_vcv_matrix(tree, ordered_names)

        C_inv = np.linalg.inv(vcv)
        ones = np.ones(n)
        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])
        Z = Y - np.outer(ones, a_hat)
        R = (Z.T @ C_inv @ Z) / (n - 1)
        eigenvalues, eigenvectors = np.linalg.eigh(R)
        idx = np.argsort(eigenvalues)[::-1]
        eigenvectors = eigenvectors[:, idx]
        scores = Z @ eigenvectors

        node_estimates, _, tree_pruned = \
            svc._reconstruct_ancestral_scores(tree, scores, ordered_names)

        name_to_idx = {name: i for i, name in enumerate(ordered_names)}
        for tip in tree_pruned.get_terminals():
            tip_scores = node_estimates[id(tip)]
            expected = scores[name_to_idx[tip.name]]
            np.testing.assert_array_almost_equal(tip_scores, expected)

    def test_root_within_tip_range(self, default_args):
        svc = PhylogeneticOrdination(default_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(MULTI_TRAITS_FILE, tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])
        vcv = svc._build_vcv_matrix(tree, ordered_names)

        C_inv = np.linalg.inv(vcv)
        ones = np.ones(n)
        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])
        Z = Y - np.outer(ones, a_hat)
        R = (Z.T @ C_inv @ Z) / (n - 1)
        eigenvalues, eigenvectors = np.linalg.eigh(R)
        idx = np.argsort(eigenvalues)[::-1]
        eigenvectors = eigenvectors[:, idx]
        scores = Z @ eigenvectors

        node_estimates, _, tree_pruned = \
            svc._reconstruct_ancestral_scores(tree, scores, ordered_names)

        root_scores = node_estimates[id(tree_pruned.root)]
        tip_pc1 = scores[:, 0]
        assert root_scores[0] >= min(tip_pc1) - 0.5
        assert root_scores[0] <= max(tip_pc1) + 0.5

    def test_node_distances_fast_path_does_not_call_distance(
        self, default_args, monkeypatch
    ):
        svc = PhylogeneticOrdination(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,C:3):0.5;"), "newick")
        ordered_names = ["A", "B", "C"]
        scores = np.array([[0.0, 0.0], [2.0, 2.0], [4.0, 4.0]])

        def fail_distance(self, *args, **kwargs):
            raise AssertionError("distance fallback should not be called")

        monkeypatch.setattr(type(tree), "distance", fail_distance)
        _, node_distances, tree_pruned = svc._reconstruct_ancestral_scores(
            tree, scores, ordered_names
        )

        assert node_distances[id(tree_pruned.root)] == pytest.approx(0.0)
        tip_distances = {
            tip.name: node_distances[id(tip)]
            for tip in tree_pruned.get_terminals()
        }
        assert tip_distances == pytest.approx({"A": 2.0, "B": 2.0, "C": 3.0})

    def test_reconstruct_ancestral_scores_uses_direct_traversal(
        self, default_args, monkeypatch
    ):
        svc = PhylogeneticOrdination(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,C:3):0.5;"), "newick")
        ordered_names = ["A", "B", "C"]
        scores = np.array([[0.0, 0.0], [2.0, 0.0], [4.0, 0.0]])

        def fail_generic_traversal(*_args, **_kwargs):
            raise AssertionError("generic tree traversal should not be used")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_generic_traversal)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_generic_traversal)
        monkeypatch.setattr(TreeMixin, "depths", fail_generic_traversal)
        monkeypatch.setattr(type(tree), "distance", fail_generic_traversal)

        node_estimates, node_distances, tree_pruned = (
            svc._reconstruct_ancestral_scores(tree, scores, ordered_names)
        )

        root = tree_pruned.root
        left_internal = root.clades[0]
        tips = {
            root.clades[0].clades[0].name: root.clades[0].clades[0],
            root.clades[0].clades[1].name: root.clades[0].clades[1],
            root.clades[1].name: root.clades[1],
        }

        np.testing.assert_allclose(node_estimates[id(left_internal)], [1.0, 0.0])
        np.testing.assert_allclose(node_estimates[id(root)], [2.0, 0.0])
        assert node_distances[id(root)] == pytest.approx(0.0)
        assert node_distances[id(tips["A"])] == pytest.approx(2.0)
        assert node_distances[id(tips["C"])] == pytest.approx(3.0)

    def test_preorder_clades_direct_handles_mixed_child_counts(
        self, default_args, monkeypatch
    ):
        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);"),
            "newick",
        )

        def fail_generic_traversal(*_args, **_kwargs):
            raise AssertionError("standard preorder traversal should be direct")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_generic_traversal)

        clades = PhylogeneticOrdination._preorder_clades_direct(tree)

        assert [clade.name for clade in clades if not clade.clades] == [
            "A", "B", "C", "D", "E", "F",
        ]
        assert clades[:3] == [
            tree.root,
            tree.root.clades[0],
            tree.root.clades[1],
        ]

    def test_draw_tree_overlay_uses_direct_preorder(
        self, default_args, monkeypatch
    ):
        pytest.importorskip("matplotlib")
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        svc = PhylogeneticOrdination(default_args)
        svc.tree_color_by = None
        tree = Phylo.read(StringIO("((A:1,B:1):1,C:3):0.5;"), "newick")
        clades = PhylogeneticOrdination._preorder_clades_direct(tree)
        node_estimates = {
            id(clade): np.array([float(index), float(index + 1)])
            for index, clade in enumerate(clades)
        }
        node_distances = {id(clade): float(index) for index, clade in enumerate(clades)}

        def fake_reconstruct(*_args, **_kwargs):
            return node_estimates, node_distances, tree

        def fail_generic_traversal(*_args, **_kwargs):
            raise AssertionError("tree overlay should use direct preorder traversal")

        monkeypatch.setattr(svc, "_reconstruct_ancestral_scores", fake_reconstruct)
        monkeypatch.setattr(TreeMixin, "find_clades", fail_generic_traversal)

        fig, ax = plt.subplots()
        try:
            svc._draw_tree_overlay(
                ax,
                fig,
                tree,
                np.zeros((3, 2)),
                ["A", "B", "C"],
                ["trait1", "trait2"],
                np.zeros((3, 2)),
            )
            assert len(ax.collections) == 1
        finally:
            plt.close(fig)


class TestParseColorBy:
    def test_column_name_match(self, default_args):
        svc = PhylogeneticOrdination(default_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(MULTI_TRAITS_FILE, tips)
        ordered_names = sorted(traits.keys())
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])

        values, categories, kind = svc._parse_color_by(
            "body_mass", trait_names, Y, ordered_names
        )
        assert kind == "continuous"
        assert len(values) == len(ordered_names)
        assert len(categories) == 0

    def test_continuous_file(self, default_args, tmp_path):
        svc = PhylogeneticOrdination(default_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(MULTI_TRAITS_FILE, tips)
        ordered_names = sorted(traits.keys())
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])

        color_file = tmp_path / "colors.tsv"
        lines = [f"{name}\t{i * 1.5}\n" for i, name in enumerate(ordered_names)]
        color_file.write_text("".join(lines))

        values, categories, kind = svc._parse_color_by(
            str(color_file), trait_names, Y, ordered_names
        )
        assert kind == "continuous"
        assert len(values) == len(ordered_names)

    def test_discrete_file(self, default_args, tmp_path):
        svc = PhylogeneticOrdination(default_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(MULTI_TRAITS_FILE, tips)
        ordered_names = sorted(traits.keys())
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])

        color_file = tmp_path / "groups.tsv"
        lines = [f"{name}\t{'A' if i % 2 == 0 else 'B'}\n"
                 for i, name in enumerate(ordered_names)]
        color_file.write_text("".join(lines))

        values, categories, kind = svc._parse_color_by(
            str(color_file), trait_names, Y, ordered_names
        )
        assert kind == "discrete"
        assert set(categories) == {"A", "B"}
        assert len(values) == len(ordered_names)

    def test_invalid_input(self, default_args):
        svc = PhylogeneticOrdination(default_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(MULTI_TRAITS_FILE, tips)
        ordered_names = sorted(traits.keys())
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])

        with pytest.raises(PhykitUserError):
            svc._parse_color_by(
                "nonexistent_column_and_not_a_file", trait_names, Y, ordered_names
            )


GENE_TREES_FILE = str(SAMPLE_FILES / "gene_trees_simple.nwk")


class TestDiscordanceVCV:
    def test_run_pca_with_gene_trees_json(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="pca",
            correction="BM",
            mode="cov",
            n_components=2,
            perplexity=None,
            n_neighbors=None,
            min_dist=0.1,
            seed=None,
            json=True,
            plot=False,
            plot_output="test_output.png",
            plot_tree=False,
            no_plot_tree=False,
            color_by=None,
            tree_color_by=None,
            gene_trees=GENE_TREES_FILE,
        )
        svc = PhylogeneticOrdination(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "eigenvalues" in data
        assert "vcv_metadata" in data
        assert data["vcv_metadata"]["n_gene_trees"] == 10

    def test_discordance_changes_pca_scores(self, capsys):
        """With discordant gene trees, PCA scores should differ."""
        base_kwargs = dict(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="pca",
            correction="BM",
            mode="cov",
            n_components=2,
            perplexity=None,
            n_neighbors=None,
            min_dist=0.1,
            seed=None,
            json=True,
            plot=False,
            plot_output="test_output.png",
            plot_tree=False,
            no_plot_tree=False,
            color_by=None,
            tree_color_by=None,
        )
        args_no_gt = Namespace(**base_kwargs)
        svc = PhylogeneticOrdination(args_no_gt)
        svc.run()
        out, _ = capsys.readouterr()
        scores_no_gt = json.loads(out)["scores"]

        args_gt = Namespace(**base_kwargs, gene_trees=GENE_TREES_FILE)
        svc = PhylogeneticOrdination(args_gt)
        svc.run()
        out, _ = capsys.readouterr()
        scores_gt = json.loads(out)["scores"]

        # At least one taxon's scores should differ
        first_taxon = sorted(scores_no_gt.keys())[0]
        pc1_no_gt = scores_no_gt[first_taxon]["PC1"]
        pc1_gt = scores_gt[first_taxon]["PC1"]
        assert pc1_gt != pytest.approx(pc1_no_gt, abs=1e-6)
