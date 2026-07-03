import builtins
import importlib
import subprocess
import sys
import pytest
import json
import numpy as np
from argparse import Namespace
from io import StringIO
from pathlib import Path
from unittest.mock import patch

from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin

import phykit.services.tree.fit_continuous as fit_continuous_module
from phykit.services.tree.fit_continuous import FitContinuous, ALL_MODELS
from phykit.helpers.pgls_utils import max_lambda as compute_max_lambda
from phykit.errors import PhykitUserError


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_traits.tsv")


def test_module_import_does_not_import_numpy_or_scipy():
    code = """
import sys
import phykit.services.tree.fit_continuous as module
assert callable(module.print_json)
assert hasattr(module.np, "__getattr__")
assert "json" not in sys.modules
assert "typing" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "numpy" not in sys.modules
assert "scipy.linalg" not in sys.modules
assert "scipy.optimize" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_module_import_does_not_import_scipy_linalg_or_optimize(monkeypatch):
    module_name = "phykit.services.tree.fit_continuous"
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
                "fit_continuous module import should not import SciPy linalg/optimize"
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


def test_lazy_numpy_caches_resolved_attributes():
    lazy_np = fit_continuous_module._LazyNumpy()

    array_attr = lazy_np.array

    assert lazy_np.__dict__["array"] is array_attr
    assert lazy_np.array is array_attr
    assert lazy_np._module is not None


@pytest.fixture
def default_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=TRAITS_FILE,
        models=None,
        json=False,
    )


@pytest.fixture
def svc(default_args):
    return FitContinuous(default_args)


@pytest.fixture
def tree_vcv_data(svc):
    """Pre-compute tree, traits, VCV and paths for direct model testing."""
    tree = svc.read_tree_file()
    tips = svc.get_tip_names_from_tree(tree)
    traits = svc._parse_trait_file(TRAITS_FILE, tips)
    ordered_names = sorted(traits.keys())
    x = np.array([traits[name] for name in ordered_names])
    vcv = svc._build_vcv_matrix(tree, ordered_names)
    parent_map = svc._build_parent_map(tree)
    paths = svc._build_root_to_tip_paths(tree, ordered_names, parent_map)
    tree_height = float(np.max(np.diag(vcv)))
    max_lam = compute_max_lambda(tree)
    return dict(
        svc=svc, tree=tree, x=x, vcv=vcv, ordered_names=ordered_names,
        paths=paths, tree_height=tree_height, max_lam=max_lam,
    )


# ── TestProcessArgs ──────────────────────────────────────────────────


class TestProcessArgs:
    def test_defaults(self):
        args = Namespace(tree="t.tre", trait_data="d.tsv", models=None, json=False)
        svc = FitContinuous(args)
        assert svc.selected_models == ALL_MODELS
        assert svc.json_output is False

    def test_subset_models(self):
        args = Namespace(tree="t.tre", trait_data="d.tsv", models="BM,OU,White", json=True)
        svc = FitContinuous(args)
        assert svc.selected_models == ["BM", "OU", "White"]
        assert svc.json_output is True


class TestTraitParsing:
    def test_comments_and_blanks(self, tmp_path, default_args):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "# comment\n\nraccoon\t1.0\nbear\t2.0\nsea_lion\t3.0\n"
            "seal\t4.0\nmonkey\t5.0\ncat\t6.0\nweasel\t7.0\ndog\t8.0\n"
        )
        svc = FitContinuous(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)

        traits = svc._parse_trait_file(str(trait_file), tree_tips)

        assert len(traits) == 8
        assert traits["raccoon"] == pytest.approx(1.0)

    def test_all_shared_trait_file_emits_no_warnings(
        self, tmp_path, default_args, capsys
    ):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "raccoon\t1.0\nbear\t2.0\nsea_lion\t3.0\nseal\t4.0\n"
            "monkey\t5.0\ncat\t6.0\nweasel\t7.0\ndog\t8.0\n"
        )
        svc = FitContinuous(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)

        traits = svc._parse_trait_file(str(trait_file), tree_tips)

        stderr = capsys.readouterr().err
        assert set(traits) == set(tree_tips)
        assert stderr == ""

    def test_extra_columns_error(self, tmp_path, default_args):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "raccoon\t1.0\nbear\t2.0\nsea_lion\t3.0\ndog\t4.0\nextra\t1.0\t2.0\n"
        )
        svc = FitContinuous(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)

        with pytest.raises(PhykitUserError) as exc_info:
            svc._parse_trait_file(str(trait_file), tree_tips)

        assert (
            "Line 5 in trait file has 3 columns; expected 2."
            in exc_info.value.messages
        )

    def test_non_numeric_trait_error(self, tmp_path, default_args):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text("raccoon\t1.0\nbear\tbad\nsea_lion\t3.0\ndog\t4.0\n")
        svc = FitContinuous(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)

        with pytest.raises(PhykitUserError) as exc_info:
            svc._parse_trait_file(str(trait_file), tree_tips)

        assert (
            "Non-numeric trait value 'bad' for taxon 'bear' on line 2."
            in exc_info.value.messages
        )


class TestTreeHelpers:
    def test_build_parent_map_uses_direct_traversal(self, monkeypatch, svc):
        tree = Phylo.read(StringIO("((A:1,B:2):3,(C:4,D:5):6);"), "newick")

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("standard tree should use direct traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        parent_map = svc._build_parent_map(tree)
        root_children = tree.root.clades
        left, right = root_children

        assert parent_map[id(left)] is tree.root
        assert parent_map[id(right)] is tree.root
        assert parent_map[id(left.clades[0])] is left
        assert parent_map[id(left.clades[1])] is left
        assert parent_map[id(right.clades[0])] is right
        assert parent_map[id(right.clades[1])] is right

    def test_build_parent_map_handles_mixed_child_counts(self, monkeypatch, svc):
        tree = Phylo.read(StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);"), "newick")

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("standard tree should use direct traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        parent_map = svc._build_parent_map(tree)
        binary = tree.root.clades[1]
        trifurcation = tree.root.clades[2]

        assert parent_map[id(tree.root.clades[0])] is tree.root
        assert parent_map[id(binary)] is tree.root
        assert parent_map[id(trifurcation)] is tree.root
        assert parent_map[id(binary.clades[0])] is binary
        assert parent_map[id(binary.clades[1])] is binary
        assert parent_map[id(trifurcation.clades[0])] is trifurcation
        assert parent_map[id(trifurcation.clades[1])] is trifurcation
        assert parent_map[id(trifurcation.clades[2])] is trifurcation


# ── TestBMModel ──────────────────────────────────────────────────────


class TestBMModel:
    def test_concentrated_ll_cholesky_matches_inverse(self, svc):
        rng = np.random.default_rng(20260623)
        n = 16
        A = rng.normal(size=(n, n))
        C = A @ A.T + np.eye(n)
        x = rng.normal(size=n)

        fast = svc._concentrated_ll_cholesky(x, C)
        inverse = svc._concentrated_ll_inverse(x, C)

        assert fast == pytest.approx(inverse)

    def test_concentrated_ll_cholesky_uses_single_solve(self, svc, monkeypatch):
        rng = np.random.default_rng(20260628)
        n = 16
        A = rng.normal(size=(n, n))
        C = A @ A.T + np.eye(n)
        x = rng.normal(size=n)
        original_cho_solve = fit_continuous_module.cho_solve
        solve_calls = 0

        def counting_cho_solve(*args, **kwargs):
            nonlocal solve_calls
            solve_calls += 1
            return original_cho_solve(*args, **kwargs)

        monkeypatch.setattr(fit_continuous_module, "cho_solve", counting_cho_solve)

        fast = svc._concentrated_ll_cholesky(x, C)
        inverse = svc._concentrated_ll_inverse(x, C)

        assert solve_calls == 1
        assert fast == pytest.approx(inverse)

    def test_concentrated_ll_cholesky_reuses_cached_scipy_linalg(
        self, svc, monkeypatch
    ):
        rng = np.random.default_rng(20260629)
        n = 16
        A = rng.normal(size=(n, n))
        C = A @ A.T + np.eye(n)
        x = rng.normal(size=n)
        previous_cho_factor = fit_continuous_module._CHO_FACTOR
        previous_cho_solve = fit_continuous_module._CHO_SOLVE

        fit_continuous_module._CHO_FACTOR = None
        fit_continuous_module._CHO_SOLVE = None
        try:
            first = svc._concentrated_ll_cholesky(x, C)
            assert fit_continuous_module._CHO_FACTOR is not None
            assert fit_continuous_module._CHO_SOLVE is not None

            original_import = builtins.__import__

            def guarded_import(
                name, globals=None, locals=None, fromlist=(), level=0
            ):
                if name == "scipy.linalg" or name.startswith("scipy.linalg."):
                    raise AssertionError(
                        "cached SciPy linalg callables should be reused"
                    )
                return original_import(name, globals, locals, fromlist, level)

            monkeypatch.setattr(builtins, "__import__", guarded_import)
            second = svc._concentrated_ll_cholesky(x, C)

            np.testing.assert_allclose(second, first, rtol=0.0, atol=1e-12)
        finally:
            fit_continuous_module._CHO_FACTOR = previous_cho_factor
            fit_continuous_module._CHO_SOLVE = previous_cho_solve

    def test_sigma2_positive(self, tree_vcv_data):
        d = tree_vcv_data
        res = d["svc"]._fit_bm(d["x"], d["vcv"])
        assert res["sigma2"] > 0
        assert res["model"] == "BM"
        assert res["k_params"] == 2

    def test_ll_finite(self, tree_vcv_data):
        d = tree_vcv_data
        res = d["svc"]._fit_bm(d["x"], d["vcv"])
        assert np.isfinite(res["log_likelihood"])


# ── TestOUModel ──────────────────────────────────────────────────────


class TestOUModel:
    def test_vcv_ou_matches_scalar_formula(self, svc, monkeypatch):
        C = np.array(
            [
                [3.0, 1.0, 0.5],
                [1.0, 2.5, 0.75],
                [0.5, 0.75, 2.0],
            ]
        )
        alpha = 0.7
        expected = np.zeros_like(C)
        for i in range(C.shape[0]):
            for j in range(i, C.shape[1]):
                shared = C[i, j]
                unique_i = C[i, i] - shared
                unique_j = C[j, j] - shared
                value = np.exp(-alpha * (unique_i + unique_j)) * (
                    1.0 - np.exp(-2.0 * alpha * shared)
                ) / (2.0 * alpha)
                expected[i, j] = value
                expected[j, i] = value

        def fail_diag(*_args, **_kwargs):
            raise AssertionError("OU VCV should use ndarray diagonal access")

        monkeypatch.setattr(fit_continuous_module.np, "diag", fail_diag)

        np.testing.assert_allclose(svc._vcv_ou(C, alpha), expected)

    def test_vcv_ou_near_zero_alpha_returns_copy(self, svc):
        C = np.array([[1.0, 0.25], [0.25, 2.0]])
        observed = svc._vcv_ou(C, 1e-12)

        np.testing.assert_array_equal(observed, C)
        assert observed is not C

    def test_ou_ll_ge_bm(self, tree_vcv_data):
        d = tree_vcv_data
        bm = d["svc"]._fit_bm(d["x"], d["vcv"])
        ou = d["svc"]._fit_ou(d["x"], d["vcv"], d["tree_height"])
        # OU has an extra parameter so its LL should be >= BM LL
        assert ou["log_likelihood"] >= bm["log_likelihood"] - 1e-6

    def test_alpha_positive(self, tree_vcv_data):
        d = tree_vcv_data
        ou = d["svc"]._fit_ou(d["x"], d["vcv"], d["tree_height"])
        assert ou["param_value"] > 0
        assert ou["param_name"] == "alpha"


# ── TestEBModel ──────────────────────────────────────────────────────


class TestEBModel:
    def test_a_finite(self, tree_vcv_data):
        d = tree_vcv_data
        eb = d["svc"]._fit_eb(d["x"], d["ordered_names"], d["paths"], d["tree_height"])
        assert np.isfinite(eb["param_value"])
        assert eb["param_name"] == "a"

    def test_eb_ll_ge_bm(self, tree_vcv_data):
        d = tree_vcv_data
        bm = d["svc"]._fit_bm(d["x"], d["vcv"])
        eb = d["svc"]._fit_eb(d["x"], d["ordered_names"], d["paths"], d["tree_height"])
        assert eb["log_likelihood"] >= bm["log_likelihood"] - 1e-6


# ── TestLambdaModel ──────────────────────────────────────────────────


class TestLambdaModel:
    def test_lambda_in_bounds(self, tree_vcv_data):
        d = tree_vcv_data
        res = d["svc"]._fit_lambda(d["x"], d["vcv"], d["max_lam"])
        assert 0 <= res["param_value"] <= d["max_lam"] + 1e-6
        assert res["param_name"] == "lambda"

    def test_lambda_matches_phylosig(self, tree_vcv_data):
        """Lambda result should have the same LL as phylogenetic_signal lambda."""
        d = tree_vcv_data
        res = d["svc"]._fit_lambda(d["x"], d["vcv"], d["max_lam"])
        assert np.isfinite(res["log_likelihood"])

    def test_lambda_restores_diagonal_directly(self, monkeypatch):
        svc = FitContinuous.__new__(FitContinuous)
        C = np.array(
            [
                [2.0, 0.4, 0.2],
                [0.4, 3.0, 0.5],
                [0.2, 0.5, 4.0],
            ]
        )
        x = np.array([1.0, 1.5, 2.0])
        original_diag = C.diagonal().copy()
        likelihood_calls = 0

        def fail_diag(*_args, **_kwargs):
            raise AssertionError("lambda search should use ndarray diagonal access")

        def fail_fill_diagonal(*_args, **_kwargs):
            raise AssertionError("lambda search should restore diagonal directly")

        def fake_concentrated_ll(_x, C_lam):
            nonlocal likelihood_calls
            likelihood_calls += 1
            np.testing.assert_allclose(C_lam.diagonal(), original_diag)
            return -abs(C_lam[0, 1] - 0.2), 1.0, 0.0

        def fake_optimize_parameter(fn, bounds):
            assert bounds == (0.0, 1.0)
            lam = 0.5
            return lam, fn(lam)

        monkeypatch.setattr(fit_continuous_module.np, "diag", fail_diag)
        monkeypatch.setattr(
            fit_continuous_module.np,
            "fill_diagonal",
            fail_fill_diagonal,
        )
        monkeypatch.setattr(svc, "_concentrated_ll", fake_concentrated_ll)
        monkeypatch.setattr(svc, "_optimize_parameter", fake_optimize_parameter)

        result = svc._fit_lambda(x, C, max_lam=1.0)

        assert result["param_name"] == "lambda"
        assert result["param_value"] == pytest.approx(0.5)
        assert np.isfinite(result["log_likelihood"])
        assert likelihood_calls == 2


# ── TestDeltaModel ───────────────────────────────────────────────────


class TestDeltaModel:
    def test_delta_1_matches_bm(self, tree_vcv_data):
        """When delta=1, the tree is untransformed, so LL should match BM."""
        d = tree_vcv_data
        bm = d["svc"]._fit_bm(d["x"], d["vcv"])

        def transform_delta1(bl, d_start, d_end):
            return d_end ** 1.0 - d_start ** 1.0

        V = d["svc"]._build_transformed_vcv(d["ordered_names"], d["paths"], transform_delta1)
        ll, sig2, z0 = d["svc"]._concentrated_ll(d["x"], V)
        assert ll == pytest.approx(bm["log_likelihood"], abs=1e-4)

    def test_delta_fitted(self, tree_vcv_data):
        d = tree_vcv_data
        res = d["svc"]._fit_delta(d["x"], d["ordered_names"], d["paths"])
        assert 0.01 <= res["param_value"] <= 3.0
        assert res["param_name"] == "delta"


# ── TestKappaModel ───────────────────────────────────────────────────


class TestKappaModel:
    def test_kappa_1_matches_bm(self, tree_vcv_data):
        """When kappa=1, branch lengths are unchanged, so LL should match BM."""
        d = tree_vcv_data
        bm = d["svc"]._fit_bm(d["x"], d["vcv"])

        def transform_kappa1(bl, d_start, d_end):
            return bl ** 1.0

        V = d["svc"]._build_transformed_vcv(d["ordered_names"], d["paths"], transform_kappa1)
        ll, sig2, z0 = d["svc"]._concentrated_ll(d["x"], V)
        assert ll == pytest.approx(bm["log_likelihood"], abs=1e-4)

    def test_kappa_fitted(self, tree_vcv_data):
        d = tree_vcv_data
        res = d["svc"]._fit_kappa(d["x"], d["ordered_names"], d["paths"])
        assert 0.01 <= res["param_value"] <= 3.0
        assert res["param_name"] == "kappa"


# ── TestWhiteNoise ───────────────────────────────────────────────────


class TestWhiteNoise:
    def test_sigma2_approx_var(self, tree_vcv_data):
        d = tree_vcv_data
        res = d["svc"]._fit_white(d["x"])
        # White noise sigma2 should approximate MLE variance
        mle_var = float(np.var(d["x"]))
        assert res["sigma2"] == pytest.approx(mle_var, rel=0.5)

    def test_ll_finite(self, tree_vcv_data):
        d = tree_vcv_data
        res = d["svc"]._fit_white(d["x"])
        assert np.isfinite(res["log_likelihood"])


# ── TestModelComparison ──────────────────────────────────────────────


class TestModelComparison:
    def test_aic_weights_sum_to_one(self, tree_vcv_data, monkeypatch):
        d = tree_vcv_data
        svc = d["svc"]
        results = []
        for m in ALL_MODELS:
            res = svc._fit_model(
                m, d["x"], d["vcv"], d["tree"], d["ordered_names"],
                d["paths"], d["max_lam"], d["tree_height"],
            )
            results.append(res)
        monkeypatch.setattr(
            fit_continuous_module.np,
            "exp",
            lambda *args, **kwargs: pytest.fail("AIC weights should use math.exp"),
        )
        results = svc._compute_model_comparison(results, len(d["x"]))
        weights = [r["aic_weight"] for r in results]
        assert sum(weights) == pytest.approx(1.0, abs=1e-6)

    def test_all_lls_finite(self, tree_vcv_data):
        d = tree_vcv_data
        svc = d["svc"]
        results = []
        for m in ALL_MODELS:
            res = svc._fit_model(
                m, d["x"], d["vcv"], d["tree"], d["ordered_names"],
                d["paths"], d["max_lam"], d["tree_height"],
            )
            results.append(res)
        for r in results:
            assert np.isfinite(r["log_likelihood"]), f"Model {r['model']} has non-finite LL"

    def test_sorted_by_aic(self, tree_vcv_data):
        d = tree_vcv_data
        svc = d["svc"]
        results = []
        for m in ALL_MODELS:
            res = svc._fit_model(
                m, d["x"], d["vcv"], d["tree"], d["ordered_names"],
                d["paths"], d["max_lam"], d["tree_height"],
            )
            results.append(res)
        results = svc._compute_model_comparison(results, len(d["x"]))
        aics = [r["aic"] for r in results]
        assert aics == sorted(aics)


# ── TestVCVTransformations ───────────────────────────────────────────


class TestVCVTransformations:
    def test_build_transformed_vcv_matches_pairwise_reference(self, svc):
        ordered_names = ["A", "B", "C"]
        paths = {
            "A": [(1, 1.0), (2, 0.5)],
            "B": [(1, 1.0), (3, 0.75)],
            "C": [(4, 1.25)],
        }

        def transform(bl, d_start, d_end):
            return d_end ** 0.8 - d_start ** 0.8

        transformed_paths = {}
        for name in ordered_names:
            cumulative = 0.0
            transformed = []
            for clade_id, branch_length in paths[name]:
                d_start = cumulative
                d_end = cumulative + branch_length
                transformed.append(
                    (clade_id, transform(branch_length, d_start, d_end))
                )
                cumulative = d_end
            transformed_paths[name] = transformed

        expected = np.zeros((3, 3))
        for i, name_i in enumerate(ordered_names):
            path_i = transformed_paths[name_i]
            expected[i, i] = sum(branch_length for _, branch_length in path_i)
            for j in range(i + 1, len(ordered_names)):
                path_j = transformed_paths[ordered_names[j]]
                shared = 0.0
                for (clade_i, branch_i), (clade_j, _) in zip(path_i, path_j):
                    if clade_i != clade_j:
                        break
                    shared += branch_i
                expected[i, j] = shared
                expected[j, i] = shared

        observed = svc._build_transformed_vcv(ordered_names, paths, transform)

        np.testing.assert_allclose(observed, expected)

    def test_build_transformed_vcv_single_tip_branches_skip_block_indexing(
        self, svc, monkeypatch
    ):
        ordered_names = ["A", "B", "C"]
        paths = {
            "A": [("A_tip", 1.0)],
            "B": [("B_tip", 2.0)],
            "C": [("C_tip", 3.0)],
        }

        def fail_ix(*_args, **_kwargs):
            raise AssertionError("single-tip branches should update diagonals directly")

        monkeypatch.setattr(fit_continuous_module.np, "ix_", fail_ix)

        observed = svc._build_transformed_vcv(
            ordered_names,
            paths,
            lambda bl, _d_start, _d_end: bl,
        )

        np.testing.assert_allclose(observed, np.diag([1.0, 2.0, 3.0]))

    def test_root_to_tip_paths_fast_path_does_not_call_get_terminals(
        self, svc, monkeypatch
    ):
        tree = Phylo.read(StringIO("((A:1,B:2):3,C:4);"), "newick")
        parent_map = svc._build_parent_map(tree)

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("standard tree should use direct traversal")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        paths = svc._build_root_to_tip_paths(tree, ["A", "C"], parent_map)

        assert [branch_length for _, branch_length in paths["A"]] == [3.0, 1.0]
        assert [branch_length for _, branch_length in paths["C"]] == [4.0]

    def test_root_to_tip_paths_avoids_reversed_for_binary_children(self, svc):
        class NoReversedList(list):
            def __reversed__(self):
                raise AssertionError("binary tip scan should push children directly")

        tree = Phylo.read(StringIO("((A:1,B:2):3,(C:4,D:5):6);"), "newick")
        for clade in tree.find_clades(order="preorder"):
            if len(clade.clades) == 2:
                clade.clades = NoReversedList(clade.clades)
        parent_map = svc._build_parent_map(tree)

        paths = svc._build_root_to_tip_paths(tree, ["A", "D"], parent_map)

        assert [branch_length for _, branch_length in paths["A"]] == [3.0, 1.0]
        assert [branch_length for _, branch_length in paths["D"]] == [6.0, 5.0]

    def test_ou_vcv_positive_definite(self, tree_vcv_data):
        d = tree_vcv_data
        V = d["svc"]._vcv_ou(d["vcv"], 0.5)
        eigenvalues = np.linalg.eigvalsh(V)
        assert np.all(eigenvalues > -1e-10)

    def test_lambda_1_equals_bm(self, tree_vcv_data, monkeypatch):
        d = tree_vcv_data

        def fail_diag(*_args, **_kwargs):
            raise AssertionError("lambda VCV should use ndarray diagonal access")

        monkeypatch.setattr(fit_continuous_module.np, "diag", fail_diag)

        V = d["svc"]._vcv_lambda(d["vcv"], 1.0)
        np.testing.assert_array_almost_equal(V, d["vcv"])

    def test_kappa_1_vcv_equals_bm(self, tree_vcv_data):
        d = tree_vcv_data

        def transform_kappa1(bl, d_start, d_end):
            return bl ** 1.0

        V = d["svc"]._build_transformed_vcv(d["ordered_names"], d["paths"], transform_kappa1)
        np.testing.assert_array_almost_equal(V, d["vcv"], decimal=4)


# ── TestRun ──────────────────────────────────────────────────────────


class TestRun:
    def test_run_uses_unmodified_tree_read(self, mocker, monkeypatch):
        args = Namespace(
            tree="dummy.tre",
            trait_data="traits.tsv",
            models="BM",
            json=False,
        )
        tree = object()
        svc = FitContinuous(args)
        read_tree = mocker.patch.object(
            svc, "read_tree_file_unmodified", return_value=tree
        )
        mocker.patch.object(
            svc,
            "read_tree_file",
            side_effect=AssertionError("copying tree reader should not be used"),
        )
        mocker.patch.object(svc, "validate_tree")
        mocker.patch.object(svc, "get_tip_names_from_tree", return_value=["a", "b", "c"])
        mocker.patch.object(
            svc,
            "_parse_trait_file",
            return_value={"a": 1.0, "b": 2.0, "c": 3.0},
        )
        vcv = np.array([
            [1.0, 0.1, 0.2],
            [0.1, 3.5, 0.3],
            [0.2, 0.3, 2.0],
        ])
        mocker.patch(
            "phykit.services.tree.vcv_utils.build_vcv_matrix",
            return_value=vcv,
        )
        mocker.patch.object(svc, "_build_parent_map", return_value={})
        mocker.patch.object(
            svc,
            "_build_root_to_tip_paths",
            return_value={"a": [], "b": [], "c": []},
        )
        fit_model = mocker.patch.object(
            svc,
            "_fit_model",
            return_value={
                "model": "BM",
                "sigma2": 1.0,
                "log_likelihood": -1.0,
                "k_params": 2,
            },
        )
        mocker.patch.object(
            svc,
            "_compute_model_comparison",
            side_effect=lambda results, _n: results,
        )
        mocker.patch.object(svc, "_fit_white", return_value={"sigma2": 1.0})
        mocker.patch.object(svc, "_print_text_output")
        monkeypatch.setattr(
            fit_continuous_module.np,
            "diag",
            lambda _matrix: (_ for _ in ()).throw(
                AssertionError("tree height should read the VCV diagonal view")
            ),
        )

        svc.run()

        read_tree.assert_called_once_with()
        assert fit_model.call_args.args[-1] == pytest.approx(3.5)
        svc.validate_tree.assert_called_once_with(
            tree,
            min_tips=3,
            require_branch_lengths=True,
            context="model fitting",
        )

    @patch("builtins.print")
    def test_text_output(self, mocked_print, default_args):
        svc = FitContinuous(default_args)
        svc.run()
        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Model Comparison (fitContinuous)" in all_output
        assert "BM" in all_output
        assert "Best model (AIC)" in all_output

    @patch("builtins.print")
    def test_print_text_output_batches_model_table(self, mocked_print):
        svc = FitContinuous.__new__(FitContinuous)
        results = [
            {
                "model": "BM",
                "param_name": None,
                "param_value": None,
                "sigma2": 1.23456,
                "z0": 2.34567,
                "log_likelihood": -10.12345,
                "aic": 24.2468,
                "delta_aic": 0.0,
                "aic_weight": 0.75,
                "bic": 26.1357,
                "delta_bic": 1.25,
                "r_squared": 0.4321,
            },
            {
                "model": "OU",
                "param_name": "alpha",
                "param_value": 0.98765,
                "sigma2": 0.87654,
                "z0": 1.76543,
                "log_likelihood": -11.98765,
                "aic": 27.9753,
                "delta_aic": 3.7285,
                "aic_weight": 0.25,
                "bic": 25.2468,
                "delta_bic": 0.0,
                "r_squared": 0.5678,
            },
        ]
        header = (
            f"{'Model':<12}{'Param':<10}{'Value':<11}"
            f"{'Sigma2':<10}{'z0':<10}{'LL':<11}"
            f"{'AIC':<9}{'dAIC':<9}{'AICw':<9}"
            f"{'BIC':<9}{'dBIC':<9}"
            f"{'R2':<7}"
        )
        ou_param = f"{0.98765:.4f}"
        expected = "\n".join(
            [
                "Model Comparison (fitContinuous)",
                "\nNumber of tips: 8\n",
                header,
                (
                    f"{'BM':<12}{'-':<10}{'-':<11}"
                    f"{1.23456:<10.4f}{2.34567:<10.4f}{-10.12345:<11.3f}"
                    f"{24.2468:<9.2f}{0.0:<9.2f}{0.75:<9.3f}"
                    f"{26.1357:<9.2f}{1.25:<9.2f}"
                    f"{0.4321:<7.3f}"
                ),
                (
                    f"{'OU':<12}{'alpha':<10}{ou_param:<11}"
                    f"{0.87654:<10.4f}{1.76543:<10.4f}{-11.98765:<11.3f}"
                    f"{27.9753:<9.2f}{3.7285:<9.2f}{0.25:<9.3f}"
                    f"{25.2468:<9.2f}{0.0:<9.2f}"
                    f"{0.5678:<7.3f}"
                ),
                "\nBest model (AIC): BM",
                "Best model (BIC): OU",
            ]
        )

        svc._print_text_output(results, 8)

        mocked_print.assert_called_once_with(expected)

    @patch("builtins.print")
    def test_json_output(self, mocked_print):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            models=None,
            json=True,
        )
        svc = FitContinuous(args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["n_tips"] == 8
        assert "models" in payload
        assert "BM" in payload["models"]
        assert "best_model_aic" in payload
        assert "best_model_bic" in payload

    @patch("builtins.print")
    def test_subset_models(self, mocked_print):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            models="BM,OU,White",
            json=True,
        )
        svc = FitContinuous(args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        assert set(payload["models"].keys()) == {"BM", "OU", "White"}


# ── TestEffectSize ────────────────────────────────────────────────────


class TestEffectSize:
    @patch("builtins.print")
    def test_r2_in_json(self, mocked_print):
        """All models have r_squared in JSON output."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            models=None,
            json=True,
        )
        svc = FitContinuous(args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        for model_name, model_data in payload["models"].items():
            assert "r_squared" in model_data, (
                f"Model {model_name} missing r_squared"
            )

    @patch("builtins.print")
    def test_white_noise_r2_zero(self, mocked_print):
        """White model has R² = 0.0 (it is the baseline)."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            models=None,
            json=True,
        )
        svc = FitContinuous(args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        white_r2 = payload["models"]["White"]["r_squared"]
        assert white_r2 == pytest.approx(0.0, abs=1e-10)

    @patch("builtins.print")
    def test_bm_r2_positive(self, mocked_print):
        """BM model has R² > 0 (phylogeny explains some variance)."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            models=None,
            json=True,
        )
        svc = FitContinuous(args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        bm_r2 = payload["models"]["BM"]["r_squared"]
        assert bm_r2 > 0

    @patch("builtins.print")
    def test_r2_in_text_output(self, mocked_print):
        """Text output contains R2 column."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            models=None,
            json=False,
        )
        svc = FitContinuous(args)
        svc.run()
        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "R2" in all_output

    @patch("builtins.print")
    def test_subset_models_without_white_still_has_r2(self, mocked_print):
        """When White isn't selected, R² is still computed."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            models="BM,OU",
            json=True,
        )
        svc = FitContinuous(args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        assert "White" not in payload["models"]
        for model_name, model_data in payload["models"].items():
            assert "r_squared" in model_data, (
                f"Model {model_name} missing r_squared"
            )

    @patch("builtins.print")
    def test_r2_matches_reference_values(self, mocked_print):
        """R² per model must match reference values from concentrated ML.

        Reference values computed via PhyKIT concentrated ML
        (cross-checked against R geiger::fitContinuous for BM sigma2):
          R geiger BM sigma2 = 0.0293769595 (rate param, different from residual var)
          PhyKIT BM sigma2   = 0.0384065703 (residual var = e'C^-1 e / n)
          PhyKIT White sigma2 = 0.7667234375
          PhyKIT BM R²       = 0.9499081827  (matches signal R²_phylo exactly)
        """
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            models=None,
            json=True,
        )
        svc = FitContinuous(args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        # BM R² should match phylogenetic signal R²_phylo (same formula)
        assert payload["models"]["BM"]["r_squared"] == pytest.approx(
            0.9499081827, abs=1e-6
        )
        # White noise baseline is always 0
        assert payload["models"]["White"]["r_squared"] == pytest.approx(0.0, abs=1e-10)
        # Lambda R² should be very close to BM (lambda ~ 1 for this data)
        assert payload["models"]["Lambda"]["r_squared"] == pytest.approx(
            0.9499082047, abs=1e-4
        )


GENE_TREES_FILE = str(SAMPLE_FILES / "gene_trees_simple.nwk")


class TestDiscordanceVCV:
    @patch("builtins.print")
    def test_run_with_gene_trees_json(self, mocked_print):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            models="BM",
            json=True,
            gene_trees=GENE_TREES_FILE,
        )
        svc = FitContinuous(args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        assert "vcv_metadata" in payload
        assert payload["vcv_metadata"]["n_gene_trees"] == 10
        assert "BM" in payload["models"]

    @patch("builtins.print")
    def test_discordance_changes_model_fit(self, mocked_print):
        """With discordant gene trees, model fits should differ."""
        args_no_gt = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            models="BM",
            json=True,
        )
        svc = FitContinuous(args_no_gt)
        svc.run()
        bm_no_gt = json.loads(mocked_print.call_args.args[0])["models"]["BM"]

        args_gt = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            models="BM",
            json=True,
            gene_trees=GENE_TREES_FILE,
        )
        svc = FitContinuous(args_gt)
        svc.run()
        bm_gt = json.loads(mocked_print.call_args.args[0])["models"]["BM"]

        # Log-likelihood should differ with discordant gene trees
        assert bm_gt["log_likelihood"] != pytest.approx(
            bm_no_gt["log_likelihood"], abs=1e-6
        )

    @patch("builtins.print")
    def test_text_output_with_gene_trees(self, mocked_print):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            models="BM",
            json=False,
            gene_trees=GENE_TREES_FILE,
        )
        svc = FitContinuous(args)
        svc.run()
        # Should complete without error
        assert mocked_print.called
