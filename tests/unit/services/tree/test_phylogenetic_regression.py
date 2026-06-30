import builtins
import importlib
import json
import subprocess
import sys

import pytest
import numpy as np
from argparse import Namespace
from pathlib import Path

from phykit.services.tree.phylogenetic_regression import PhylogeneticRegression
import phykit.services.tree.phylogenetic_regression as phylogenetic_regression_module
from phykit.helpers.pgls_utils import (
    max_lambda as compute_max_lambda,
    estimate_lambda,
    fit_gls,
)
from phykit.helpers.trait_parsing import parse_multi_trait_file
from phykit.errors import PhykitUserError


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
MULTI_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_multi_traits.tsv")


def test_module_import_does_not_import_scipy_linalg(monkeypatch):
    module_name = "phykit.services.tree.phylogenetic_regression"
    previous = sys.modules.pop(module_name, None)
    original_import = builtins.__import__

    def guarded_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "scipy.linalg" or name.startswith("scipy.linalg."):
            raise AssertionError(
                "phylogenetic_regression module import should not import SciPy linalg"
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
import phykit.services.tree.phylogenetic_regression as module
assert hasattr(module.np, "__getattr__")
assert module._CHO_FACTOR is None
assert module._CHO_SOLVE is None
assert callable(module.print_json)
assert callable(module.parse_multi_trait_file)
assert "typing" not in sys.modules
assert "numpy" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.trait_parsing" not in sys.modules
assert "scipy.linalg" not in sys.modules
assert "phykit.helpers.pgls_utils" not in sys.modules
    """
    subprocess.run([sys.executable, "-c", code], check=True)


def test_repeated_fit_model_caches_scipy_linalg_imports(monkeypatch):
    previous_cho_factor = phylogenetic_regression_module._CHO_FACTOR
    previous_cho_solve = phylogenetic_regression_module._CHO_SOLVE
    phylogenetic_regression_module._CHO_FACTOR = None
    phylogenetic_regression_module._CHO_SOLVE = None
    original_import = builtins.__import__
    scipy_linalg_imports = 0

    def counting_import(name, globals=None, locals=None, fromlist=(), level=0):
        nonlocal scipy_linalg_imports
        if name == "scipy.linalg":
            scipy_linalg_imports += 1
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", counting_import)
    try:
        rng = np.random.default_rng(20260628)
        svc = PhylogeneticRegression.__new__(PhylogeneticRegression)
        n = 8
        k = 2
        a_matrix = rng.normal(size=(n, n))
        vcv = a_matrix @ a_matrix.T + np.eye(n)
        X = np.column_stack([np.ones(n), rng.normal(size=(n, k))])
        y = X @ np.array([0.25, 1.1, -0.4]) + rng.normal(size=n) * 0.1

        svc._fit_model(y, X, vcv, k, n)
        first_call_imports = scipy_linalg_imports
        svc._fit_model(y, X, vcv, k, n)
    finally:
        phylogenetic_regression_module._CHO_FACTOR = previous_cho_factor
        phylogenetic_regression_module._CHO_SOLVE = previous_cho_solve

    assert first_call_imports > 0
    assert scipy_linalg_imports == first_call_imports


def test_probability_helpers_match_expected_values():
    p_values = phylogenetic_regression_module._t_two_tailed_p_values(
        np.array([0.0, 2.2281388519649385, -2.2281388519649385]),
        10,
    )
    assert p_values == pytest.approx([1.0, 0.05, 0.05])
    assert phylogenetic_regression_module._f_sf(
        2.2281388519649385 ** 2, 1, 10
    ) == pytest.approx(0.05)
    assert phylogenetic_regression_module._f_sf(4.9646027437307145, 2, 10) == pytest.approx(0.031809004261296604)


def test_probability_helpers_do_not_import_scipy_stats(monkeypatch):
    original_import = __import__

    def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "scipy.stats" or name.startswith("scipy.stats."):
            raise AssertionError("phylogenetic regression p-values should not import scipy.stats")
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr("builtins.__import__", fake_import)

    p_values = phylogenetic_regression_module._t_two_tailed_p_values(
        np.array([0.0, 2.2281388519649385]),
        10,
    )
    assert p_values == pytest.approx([1.0, 0.05])
    assert phylogenetic_regression_module._f_sf(4.9646027437307145, 2, 10) == pytest.approx(0.031809004261296604)


def test_single_predictor_f_p_value_does_not_import_fdtrc(monkeypatch):
    previous_stdtr = phylogenetic_regression_module._STDTR
    previous_fdtrc = phylogenetic_regression_module._FDTRC
    phylogenetic_regression_module._STDTR = lambda df, x: 0.025
    phylogenetic_regression_module._FDTRC = None
    original_import = builtins.__import__

    def fail_special_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "scipy.special":
            raise AssertionError("dfn=1 F p-values should not import fdtrc")
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", fail_special_import)
    try:
        assert phylogenetic_regression_module._f_sf(4.0, 1, 10) == pytest.approx(0.05)
    finally:
        phylogenetic_regression_module._STDTR = previous_stdtr
        phylogenetic_regression_module._FDTRC = previous_fdtrc


def test_probability_helpers_cache_scipy_special_imports(monkeypatch):
    previous_stdtr = phylogenetic_regression_module._STDTR
    previous_fdtrc = phylogenetic_regression_module._FDTRC
    phylogenetic_regression_module._STDTR = None
    phylogenetic_regression_module._FDTRC = None
    original_import = builtins.__import__
    special_imports = 0

    def counting_import(name, globals=None, locals=None, fromlist=(), level=0):
        nonlocal special_imports
        if name == "scipy.special":
            special_imports += 1
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", counting_import)
    try:
        phylogenetic_regression_module._t_two_tailed_p_values(
            np.array([0.0, 2.2281388519649385]),
            10,
        )
        phylogenetic_regression_module._t_two_tailed_p_values(
            np.array([1.0, -1.0]),
            10,
        )
        phylogenetic_regression_module._f_sf(4.9646027437307145, 2, 10)
        phylogenetic_regression_module._f_sf(2.978237016082323, 2, 10)
    finally:
        phylogenetic_regression_module._STDTR = previous_stdtr
        phylogenetic_regression_module._FDTRC = previous_fdtrc

    assert special_imports == 2


def test_format_result_maps_parallel_arrays_by_name():
    svc = PhylogeneticRegression.__new__(PhylogeneticRegression)

    result = svc._format_result(
        coef_names=["(Intercept)", "body_mass", "longevity"],
        beta_hat=np.array([0.25, 1.5, -0.75]),
        se=np.array([0.1, 0.2, 0.3]),
        t_stats=np.array([2.5, 7.5, -2.5]),
        p_values=np.array([0.03, 0.001, 0.04]),
        sigma2=0.64,
        df_resid=7,
        r_squared=0.8,
        adj_r_squared=0.7,
        f_stat=12.3,
        f_p_value=0.004,
        ll=-3.5,
        aic=15.0,
        formula="brain_size ~ body_mass + longevity",
        k=2,
        n=3,
        residuals=np.array([0.11, -0.22, 0.33]),
        fitted=np.array([1.1, 2.2, 3.3]),
        ordered_names=["taxon_a", "taxon_b", "taxon_c"],
        lambda_val=0.5,
        r2_total=0.9,
        r2_phylo=0.1,
    )

    assert result["coefficients"]["body_mass"] == {
        "estimate": 1.5,
        "std_error": 0.2,
        "t_value": 7.5,
        "p_value": 0.001,
    }
    assert result["residuals"] == {
        "taxon_a": 0.11,
        "taxon_b": -0.22,
        "taxon_c": 0.33,
    }
    assert result["fitted_values"] == {
        "taxon_a": 1.1,
        "taxon_b": 2.2,
        "taxon_c": 3.3,
    }
    assert result["lambda"] == 0.5


@pytest.fixture
def default_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=MULTI_TRAITS_FILE,
        response="brain_size",
        predictors=["body_mass"],
        method="BM",
        json=False,
    )


@pytest.fixture
def lambda_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=MULTI_TRAITS_FILE,
        response="brain_size",
        predictors=["body_mass"],
        method="lambda",
        json=False,
    )


@pytest.fixture
def multi_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=MULTI_TRAITS_FILE,
        response="brain_size",
        predictors=["body_mass", "longevity"],
        method="BM",
        json=False,
    )


@pytest.fixture
def json_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=MULTI_TRAITS_FILE,
        response="brain_size",
        predictors=["body_mass"],
        method="BM",
        json=True,
    )


class TestProcessArgs:
    def test_defaults(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            response="y",
            predictors=["x1"],
        )
        svc = PhylogeneticRegression.__new__(PhylogeneticRegression)
        parsed = svc.process_args(args)
        assert parsed["tree_file_path"] == "t.tre"
        assert parsed["trait_data_path"] == "d.tsv"
        assert parsed["response"] == "y"
        assert parsed["predictors"] == ["x1"]
        assert parsed["method"] == "BM"
        assert parsed["json_output"] is False

    def test_overrides(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            response="y",
            predictors=["x1", "x2"],
            method="lambda",
            json=True,
        )
        svc = PhylogeneticRegression.__new__(PhylogeneticRegression)
        parsed = svc.process_args(args)
        assert parsed["method"] == "lambda"
        assert parsed["json_output"] is True
        assert parsed["predictors"] == ["x1", "x2"]


class TestValidation:
    def test_missing_response_column(self):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            response="nonexistent",
            predictors=["body_mass"],
            method="BM",
            json=False,
        )
        with pytest.raises(SystemExit):
            svc = PhylogeneticRegression(args)
            svc.run()

    def test_missing_predictor_column(self):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            response="brain_size",
            predictors=["nonexistent"],
            method="BM",
            json=False,
        )
        with pytest.raises(SystemExit):
            svc = PhylogeneticRegression(args)
            svc.run()

    def test_response_in_predictors(self):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            response="brain_size",
            predictors=["brain_size"],
            method="BM",
            json=False,
        )
        with pytest.raises(SystemExit):
            svc = PhylogeneticRegression(args)
            svc.run()


class TestPGLSBM:
    """Test BM regression against R reference values.

    Reference values computed in R 4.4.0 using manual GLS with ape::vcv()
    (raw VCV, matching caper::pgls behavior):

        library(ape)
        tree <- read.tree("tree_simple.tre")
        traits <- read.delim("tree_simple_multi_traits.tsv", row.names=1)
        ord <- sort(rownames(traits)); traits_s <- traits[ord, ]
        C <- vcv(tree)[ord, ord]; C_inv <- solve(C)
        y <- traits_s$brain_size; X <- cbind(1, traits_s$body_mass)
        beta <- solve(t(X)%*%C_inv%*%X) %*% t(X)%*%C_inv%*%y
        # Intercept=0.9972138  body_mass=0.7085628
    """

    def test_coefficients_match_r(self, default_args):
        """Coefficients must match R manual GLS with raw VCV."""
        svc = PhylogeneticRegression(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(svc.trait_data_path, tree_tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        resp_idx = trait_names.index("brain_size")
        pred_idx = trait_names.index("body_mass")

        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X = np.column_stack([
            np.ones(n),
            np.array([traits[name][pred_idx] for name in ordered_names]),
        ])
        vcv = svc._build_vcv_matrix(tree, ordered_names)
        C_inv = np.linalg.inv(vcv)

        beta_hat, residuals, sigma2, var_beta = fit_gls(y, X, C_inv)
        se = np.sqrt(np.diag(var_beta))

        # R reference: Intercept=0.9972138 SE=0.0871177 t=11.4467407
        assert beta_hat[0] == pytest.approx(0.9972139, abs=1e-4)
        assert se[0] == pytest.approx(0.0871177, abs=1e-4)

        # R reference: body_mass=0.7085628 SE=0.0450911 t=15.7140257
        assert beta_hat[1] == pytest.approx(0.7085628, abs=1e-4)
        assert se[1] == pytest.approx(0.0450911, abs=1e-4)

    def test_model_stats_match_r(self, default_args):
        """R², F-stat, log-lik, AIC must match R manual GLS."""
        svc = PhylogeneticRegression(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(svc.trait_data_path, tree_tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        resp_idx = trait_names.index("brain_size")
        pred_idx = trait_names.index("body_mass")

        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X = np.column_stack([
            np.ones(n),
            np.array([traits[name][pred_idx] for name in ordered_names]),
        ])
        vcv = svc._build_vcv_matrix(tree, ordered_names)
        C_inv = np.linalg.inv(vcv)

        beta_hat, residuals, sigma2, var_beta = fit_gls(y, X, C_inv)
        fitted = X @ beta_hat

        r2, adj_r2, f_stat, f_p, r2_total, r2_phylo = svc._compute_model_stats(
            y, fitted, residuals, C_inv, k=1, n=n
        )

        # R reference
        assert r2 == pytest.approx(0.9762781, abs=1e-4)
        assert adj_r2 == pytest.approx(0.9723244, abs=1e-4)
        assert f_stat == pytest.approx(246.9306, abs=0.01)
        assert f_p == pytest.approx(4.2092e-06, rel=0.01)

        # Residual SE (REML): 0.0249942
        rse = np.sqrt(sigma2)
        assert rse == pytest.approx(0.0249942, abs=1e-4)

    def test_log_likelihood_and_aic_match_r(self, default_args):
        """ML log-likelihood and AIC must match R."""
        svc = PhylogeneticRegression(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(svc.trait_data_path, tree_tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        resp_idx = trait_names.index("brain_size")
        pred_idx = trait_names.index("body_mass")

        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X = np.column_stack([
            np.ones(n),
            np.array([traits[name][pred_idx] for name in ordered_names]),
        ])
        vcv = svc._build_vcv_matrix(tree, ordered_names)
        C_inv = np.linalg.inv(vcv)
        beta_hat, residuals, _, _ = fit_gls(y, X, C_inv)

        # Compute ML log-likelihood (same formula as in run())
        sign, logdet_C = np.linalg.slogdet(vcv)
        sigma2_ml = float(residuals @ C_inv @ residuals) / n
        ll = -0.5 * (n * np.log(2 * np.pi) + n * np.log(sigma2_ml) + logdet_C + n)
        aic = -2.0 * ll + 2.0 * (1 + 2)

        # R reference: ML log-likelihood=6.0558469, AIC=-6.1116939
        assert ll == pytest.approx(6.0558, abs=0.001)
        assert aic == pytest.approx(-6.1117, abs=0.001)

    def test_fit_model_cholesky_matches_inverse(self):
        rng = np.random.default_rng(20260623)
        svc = PhylogeneticRegression.__new__(PhylogeneticRegression)
        n = 16
        k = 2
        A = rng.normal(size=(n, n))
        vcv = A @ A.T + np.eye(n)
        X = np.column_stack([np.ones(n), rng.normal(size=(n, k))])
        beta = np.array([0.25, 1.1, -0.4])
        y = X @ beta + rng.normal(size=n) * 0.1

        fast = svc._fit_model(y, X, vcv, k, n)
        inverse = svc._fit_model_inverse(y, X, vcv, k, n)

        for fast_value, inverse_value in zip(fast, inverse):
            assert fast_value == pytest.approx(inverse_value)

    def test_fit_model_uses_single_cholesky_solve(self, monkeypatch):
        rng = np.random.default_rng(20260628)
        svc = PhylogeneticRegression.__new__(PhylogeneticRegression)
        n = 16
        k = 2
        a_matrix = rng.normal(size=(n, n))
        vcv = a_matrix @ a_matrix.T + np.eye(n)
        X = np.column_stack([np.ones(n), rng.normal(size=(n, k))])
        beta = np.array([0.25, 1.1, -0.4])
        y = X @ beta + rng.normal(size=n) * 0.1
        original_cho_solve = phylogenetic_regression_module.cho_solve
        solve_calls = 0

        def counting_cho_solve(*args, **kwargs):
            nonlocal solve_calls
            solve_calls += 1
            return original_cho_solve(*args, **kwargs)

        monkeypatch.setattr(
            phylogenetic_regression_module,
            "cho_solve",
            counting_cho_solve,
        )

        fast = svc._fit_model(y, X, vcv, k, n)
        inverse = svc._fit_model_inverse(y, X, vcv, k, n)

        assert solve_calls == 1
        for fast_value, inverse_value in zip(fast, inverse):
            assert fast_value == pytest.approx(inverse_value)

    def test_fit_model_cholesky_avoids_normal_matrix_inverse(self, monkeypatch):
        rng = np.random.default_rng(20260628)
        svc = PhylogeneticRegression.__new__(PhylogeneticRegression)
        n = 16
        k = 2
        a_matrix = rng.normal(size=(n, n))
        vcv = a_matrix @ a_matrix.T + np.eye(n)
        X = np.column_stack([np.ones(n), rng.normal(size=(n, k))])
        beta = np.array([0.25, 1.1, -0.4])
        y = X @ beta + rng.normal(size=n) * 0.1

        def fail_inv(*_args, **_kwargs):
            raise AssertionError("fast Cholesky path should use solve")

        monkeypatch.setattr(phylogenetic_regression_module.np.linalg, "inv", fail_inv)

        beta_hat, residuals, sigma2, var_beta, *_ = svc._fit_model(
            y, X, vcv, k, n
        )

        assert beta_hat.shape == (k + 1,)
        assert residuals.shape == (n,)
        assert sigma2 > 0
        assert var_beta.shape == (k + 1, k + 1)

    def test_t_stats_and_pvalues(self, default_args):
        svc = PhylogeneticRegression(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(svc.trait_data_path, tree_tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        resp_idx = trait_names.index("brain_size")
        pred_idx = trait_names.index("body_mass")

        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X = np.column_stack([
            np.ones(n),
            np.array([traits[name][pred_idx] for name in ordered_names]),
        ])
        vcv = svc._build_vcv_matrix(tree, ordered_names)
        C_inv = np.linalg.inv(vcv)

        beta_hat, residuals, sigma2, var_beta = fit_gls(y, X, C_inv)
        se = np.sqrt(np.diag(var_beta))
        t_stats = beta_hat / se

        from scipy.stats import t as t_distribution
        df_resid = n - 2
        p_values = 2.0 * t_distribution.sf(np.abs(t_stats), df=df_resid)

        # R reference: t-values 11.4467407, 15.7140257
        assert t_stats[0] == pytest.approx(11.4467, abs=0.01)
        assert t_stats[1] == pytest.approx(15.7140, abs=0.01)

        # R reference: p-values 0.0000266772, 0.0000042092
        assert p_values[0] == pytest.approx(2.6677e-05, rel=0.01)
        assert p_values[1] == pytest.approx(4.2092e-06, rel=0.01)


class TestPGLSLambda:
    def test_lambda_estimation(self, lambda_args):
        svc = PhylogeneticRegression(lambda_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(svc.trait_data_path, tree_tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        resp_idx = trait_names.index("brain_size")
        pred_idx = trait_names.index("body_mass")

        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X = np.column_stack([
            np.ones(n),
            np.array([traits[name][pred_idx] for name in ordered_names]),
        ])
        vcv = svc._build_vcv_matrix(tree, ordered_names)
        max_lam = compute_max_lambda(tree)

        lambda_val, ll = estimate_lambda(y, X, vcv, max_lam)

        # Lambda should be between 0 and max_lambda
        assert 0 <= lambda_val <= max_lam
        # Log-likelihood should be a finite number
        assert np.isfinite(ll)

    def test_lambda_coefficients(self, lambda_args):
        svc = PhylogeneticRegression(lambda_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(svc.trait_data_path, tree_tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        resp_idx = trait_names.index("brain_size")
        pred_idx = trait_names.index("body_mass")

        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X = np.column_stack([
            np.ones(n),
            np.array([traits[name][pred_idx] for name in ordered_names]),
        ])
        vcv = svc._build_vcv_matrix(tree, ordered_names)
        max_lam = compute_max_lambda(tree)

        lambda_val, _ = estimate_lambda(y, X, vcv, max_lam)

        # Transform VCV
        diag_vals = np.diag(vcv).copy()
        vcv_t = vcv * lambda_val
        np.fill_diagonal(vcv_t, diag_vals)

        C_inv = np.linalg.inv(vcv_t)
        beta_hat, _, _, _ = fit_gls(y, X, C_inv)

        # Coefficients should still be reasonable
        assert 0.5 < beta_hat[0] < 1.5
        assert 0.4 < beta_hat[1] < 1.0


class TestMultipleRegression:
    """Test multiple regression against R reference values.

    R reference (manual GLS with raw VCV):
        Intercept=0.6544307 SE=0.3311317 t=1.9763454 p=0.1050679
        body_mass=0.6193379 SE=0.0943992 t=6.5608407 p=0.0012333
        longevity=0.1628710 SE=0.1519291 t=1.0720198 p=0.3327044
        Residual SE (REML)=0.0246890
        ML log-likelihood=6.8834005  AIC=-5.7668009
    """

    def test_multiple_predictors_match_r(self, multi_args):
        svc = PhylogeneticRegression(multi_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(svc.trait_data_path, tree_tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        resp_idx = trait_names.index("brain_size")
        pred_indices = [trait_names.index("body_mass"), trait_names.index("longevity")]

        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X = np.column_stack([
            np.ones(n),
            np.array([[traits[name][j] for j in pred_indices] for name in ordered_names]),
        ])
        vcv = svc._build_vcv_matrix(tree, ordered_names)
        C_inv = np.linalg.inv(vcv)

        beta_hat, residuals, sigma2, var_beta = fit_gls(y, X, C_inv)
        se = np.sqrt(np.diag(var_beta))

        # Should have 3 coefficients: intercept + 2 predictors
        assert len(beta_hat) == 3
        assert var_beta.shape == (3, 3)

        # R reference values
        assert beta_hat[0] == pytest.approx(0.6544307, abs=1e-4)
        assert beta_hat[1] == pytest.approx(0.6193379, abs=1e-4)
        assert beta_hat[2] == pytest.approx(0.1628710, abs=1e-4)

        assert se[0] == pytest.approx(0.3311317, abs=1e-4)
        assert se[1] == pytest.approx(0.0943992, abs=1e-4)
        assert se[2] == pytest.approx(0.1519291, abs=1e-4)

        rse = np.sqrt(sigma2)
        assert rse == pytest.approx(0.024689, abs=1e-4)


class TestRun:
    def test_run_uses_unmodified_tree_read(self, mocker):
        args = Namespace(
            tree="/some/path/to/file.tre",
            trait_data="/some/path/to/traits.tsv",
            response="y",
            predictors=["x"],
            method="BM",
            json=False,
            gene_trees=None,
        )
        svc = PhylogeneticRegression(args)
        tree = object()
        mocked_read = mocker.patch.object(
            PhylogeneticRegression,
            "read_tree_file_unmodified",
            return_value=tree,
        )
        mocker.patch.object(
            PhylogeneticRegression,
            "read_tree_file",
            side_effect=AssertionError("run should not copy cached trees"),
        )
        mocked_validate = mocker.patch.object(
            PhylogeneticRegression, "validate_tree"
        )
        mocker.patch.object(
            PhylogeneticRegression,
            "get_tip_names_from_tree",
            return_value=["a", "b", "c"],
        )
        mocker.patch(
            "phykit.services.tree.phylogenetic_regression.parse_multi_trait_file",
            return_value=(
                ["y", "x"],
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
            PhylogeneticRegression,
            "_fit_model",
            return_value=(
                np.array([1.0, 0.5]),
                np.zeros(3),
                1.0,
                np.eye(2),
                0.5,
                0.4,
                2.0,
                0.1,
                0.6,
                0.1,
                -5.0,
                14.0,
            ),
        )
        mocker.patch(
            "phykit.services.tree.phylogenetic_regression._t_two_tailed_p_values",
            return_value=np.array([0.2, 0.3]),
        )
        mocker.patch.object(PhylogeneticRegression, "_print_text_output")

        svc.run()

        mocked_read.assert_called_once_with()
        mocked_validate.assert_called_once_with(
            tree,
            min_tips=3,
            require_branch_lengths=True,
            context="phylogenetic regression",
        )
        mocked_vcv.assert_called_once_with(tree, ["a", "b", "c"])

    def test_run_resolves_trait_columns_with_first_index_map(self, mocker):
        args = Namespace(
            tree="/some/path/to/file.tre",
            trait_data="/some/path/to/traits.tsv",
            response="y",
            predictors=["x"],
            method="BM",
            json=False,
            gene_trees=None,
        )
        svc = PhylogeneticRegression(args)
        tree = object()
        mocker.patch.object(svc, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(svc, "validate_tree")
        mocker.patch.object(svc, "get_tip_names_from_tree", return_value=["a", "b", "c"])
        mocker.patch(
            "phykit.services.tree.phylogenetic_regression.parse_multi_trait_file",
            return_value=(
                ["y", "x", "y", "x"],
                {
                    "a": [0.0, 1.0, 9.0, 9.0],
                    "b": [1.0, 2.0, 9.0, 9.0],
                    "c": [2.0, 3.0, 9.0, 9.0],
                },
            ),
        )
        mocker.patch(
            "phykit.services.tree.vcv_utils.build_vcv_matrix",
            return_value=np.eye(3),
        )

        def arrays(traits, ordered_names, response_index, predictor_indices):
            assert response_index == 0
            assert predictor_indices == [1]
            return np.array([0.0, 1.0, 2.0]), np.ones((3, 2))

        mocker.patch(
            "phykit.services.tree.phylogenetic_regression.response_predictor_arrays",
            side_effect=arrays,
        )
        mocker.patch.object(
            PhylogeneticRegression,
            "_fit_model",
            return_value=(
                np.array([1.0, 0.5]),
                np.zeros(3),
                1.0,
                np.eye(2),
                0.5,
                0.4,
                2.0,
                0.1,
                0.6,
                0.1,
                -5.0,
                14.0,
            ),
        )
        mocker.patch(
            "phykit.services.tree.phylogenetic_regression._t_two_tailed_p_values",
            return_value=np.array([0.2, 0.3]),
        )
        mocker.patch.object(PhylogeneticRegression, "_print_text_output")

        svc.run()

    def test_text_output(self, default_args, capsys):
        svc = PhylogeneticRegression(default_args)
        svc.run()
        captured = capsys.readouterr()
        assert "Phylogenetic Generalized Least Squares (PGLS)" in captured.out
        assert "Formula: brain_size ~ body_mass" in captured.out
        assert "(Intercept)" in captured.out
        assert "body_mass" in captured.out
        assert "R-squared" in captured.out
        assert "F-statistic" in captured.out
        assert "Log-likelihood" in captured.out
        assert "AIC" in captured.out

    def test_print_text_output_batches_coefficients(self, default_args, mocker):
        svc = PhylogeneticRegression(default_args)
        printed = mocker.patch("builtins.print")

        svc._print_text_output(
            coef_names=["(Intercept)", "body_mass", "longevity"],
            beta_hat=np.array([0.6544307, 0.6193379, 0.162871]),
            se=np.array([0.3311317, 0.0943992, 0.1519291]),
            t_stats=np.array([1.9763454, 6.5608407, 1.0720198]),
            p_values=np.array([0.1050679, 0.0012333, 0.3327044]),
            sigma2=0.024689,
            df_resid=5,
            r_squared=0.9,
            adj_r_squared=0.86,
            f_stat=12.34,
            f_p_value=0.00456,
            ll=6.8834005,
            aic=-5.7668009,
            formula="brain_size ~ body_mass + longevity",
            k=2,
            n=8,
            lambda_val=0.75,
            r2_total=0.92,
            r2_phylo=0.12,
        )

        printed.assert_called_once_with(
            "Phylogenetic Generalized Least Squares (PGLS)\n"
            "\n"
            "Formula: brain_size ~ body_mass + longevity\n"
            "\n"
            "Estimated lambda: 0.7500\n"
            "\n"
            "Coefficients:\n"
            "                        Estimate   Std.Error     t-value     p-value\n"
            "(Intercept)               0.6544      0.3311      1.9763    0.105068     \n"
            "body_mass                 0.6193      0.0944      6.5608    0.001233    **\n"
            "longevity                 0.1629      0.1519      1.0720    0.332704     \n"
            "---\n"
            "Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1\n"
            "\n"
            "Residual standard error: 0.1571 on 5 degrees of freedom\n"
            "Multiple R-squared: 0.9000    Adjusted R-squared: 0.8600\n"
            "R-squared (total):   0.9200   (phylo + predictor)\n"
            "R-squared (phylo):   0.1200   (phylogeny contribution)\n"
            "F-statistic: 12.34 on 2 and 5 DF    p-value: 0.004560\n"
            "Log-likelihood: 6.8834    AIC: -5.7668"
        )

    def test_json_output(self, json_args, capsys):
        svc = PhylogeneticRegression(json_args)
        svc.run()
        captured = capsys.readouterr()
        payload = json.loads(captured.out)
        assert "formula" in payload
        assert "coefficients" in payload
        assert "(Intercept)" in payload["coefficients"]
        assert "body_mass" in payload["coefficients"]
        assert "r_squared" in payload
        assert "adj_r_squared" in payload
        assert "f_statistic" in payload
        assert "log_likelihood" in payload
        assert "aic" in payload
        assert "residuals" in payload
        assert "fitted_values" in payload

    def test_format_result_uses_bulk_residual_and_fitted_conversion(self):
        class ToListOnly:
            def __init__(self, values):
                self.values = values

            def tolist(self):
                return self.values

            def __iter__(self):
                raise AssertionError("result vectors should be converted in bulk")

        svc = PhylogeneticRegression.__new__(PhylogeneticRegression)

        result = svc._format_result(
            coef_names=["(Intercept)", "body_mass"],
            beta_hat=[1.0, 2.0],
            se=[0.5, 0.25],
            t_stats=[2.0, 8.0],
            p_values=[0.1, 0.01],
            sigma2=1.0,
            df_resid=3,
            r_squared=0.2,
            adj_r_squared=0.1,
            f_stat=4.0,
            f_p_value=0.05,
            ll=-10.0,
            aic=24.0,
            formula="brain_size ~ body_mass",
            k=1,
            n=5,
            residuals=ToListOnly([0.1, -0.2]),
            fitted=ToListOnly([1.1, 1.2]),
            ordered_names=["A", "B"],
            lambda_val=None,
            r2_total=0.3,
            r2_phylo=0.1,
        )

        assert result["residuals"] == {"A": 0.1, "B": -0.2}
        assert result["fitted_values"] == {"A": 1.1, "B": 1.2}

    def test_json_lambda_output(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            response="brain_size",
            predictors=["body_mass"],
            method="lambda",
            json=True,
        )
        svc = PhylogeneticRegression(args)
        svc.run()
        captured = capsys.readouterr()
        payload = json.loads(captured.out)
        assert "lambda" in payload

    def test_text_lambda_output(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            response="brain_size",
            predictors=["body_mass"],
            method="lambda",
            json=False,
        )
        svc = PhylogeneticRegression(args)
        svc.run()
        captured = capsys.readouterr()
        assert "Estimated lambda" in captured.out

    def test_signif_codes_displayed(self, default_args, capsys):
        svc = PhylogeneticRegression(default_args)
        svc.run()
        captured = capsys.readouterr()
        assert "Signif. codes:" in captured.out
        assert "***" in captured.out

    def test_multiple_predictors_output(self, multi_args, capsys):
        svc = PhylogeneticRegression(multi_args)
        svc.run()
        captured = capsys.readouterr()
        assert "body_mass" in captured.out
        assert "longevity" in captured.out
        assert "brain_size ~ body_mass + longevity" in captured.out


GENE_TREES_FILE = str(SAMPLE_FILES / "gene_trees_simple.nwk")


class TestDiscordanceVCV:
    def test_run_with_gene_trees_json(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            response="brain_size",
            predictors=["body_mass"],
            method="BM",
            json=True,
            gene_trees=GENE_TREES_FILE,
        )
        svc = PhylogeneticRegression(args)
        svc.run()
        captured = capsys.readouterr()
        payload = json.loads(captured.out)
        assert "vcv_metadata" in payload
        assert payload["vcv_metadata"]["n_gene_trees"] == 10
        assert "coefficients" in payload

    def test_discordance_changes_coefficients(self, capsys):
        """With discordant gene trees, PGLS coefficients should differ."""
        args_no_gt = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            response="brain_size",
            predictors=["body_mass"],
            method="BM",
            json=True,
        )
        svc = PhylogeneticRegression(args_no_gt)
        svc.run()
        out, _ = capsys.readouterr()
        coefs_no_gt = json.loads(out)["coefficients"]

        args_gt = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            response="brain_size",
            predictors=["body_mass"],
            method="BM",
            json=True,
            gene_trees=GENE_TREES_FILE,
        )
        svc = PhylogeneticRegression(args_gt)
        svc.run()
        out, _ = capsys.readouterr()
        coefs_gt = json.loads(out)["coefficients"]

        # Coefficients should differ (not exactly equal)
        assert coefs_gt["body_mass"]["estimate"] != pytest.approx(
            coefs_no_gt["body_mass"]["estimate"], abs=1e-6
        )

    def test_text_output_with_gene_trees(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            response="brain_size",
            predictors=["body_mass"],
            method="BM",
            json=False,
            gene_trees=GENE_TREES_FILE,
        )
        svc = PhylogeneticRegression(args)
        svc.run()
        captured = capsys.readouterr()
        assert "PGLS" in captured.out
        assert "body_mass" in captured.out


class TestEffectSize:
    """Tests for the three-way R² variance decomposition."""

    def test_r2_phylo_in_json(self, capsys):
        """JSON output contains r_squared_phylo, r_squared_total, and r_squared."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            response="brain_size",
            predictors=["body_mass"],
            method="BM",
            json=True,
        )
        svc = PhylogeneticRegression(args)
        svc.run()
        captured = capsys.readouterr()
        payload = json.loads(captured.out)
        assert "r_squared" in payload
        assert "r_squared_phylo" in payload
        assert "r_squared_total" in payload

    def test_decomposition_sums(self, capsys):
        """R²_phylo + R²_pred should approximately equal R²_total."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            response="brain_size",
            predictors=["body_mass"],
            method="BM",
            json=True,
        )
        svc = PhylogeneticRegression(args)
        svc.run()
        captured = capsys.readouterr()
        payload = json.loads(captured.out)
        r2_pred = payload["r_squared"]
        r2_total = payload["r_squared_total"]
        r2_phylo = payload["r_squared_phylo"]
        assert r2_phylo + r2_pred == pytest.approx(r2_total, abs=1e-10)

    def test_r2_values_finite(self, capsys):
        """Both new R² values should be finite."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            response="brain_size",
            predictors=["body_mass"],
            method="BM",
            json=True,
        )
        svc = PhylogeneticRegression(args)
        svc.run()
        captured = capsys.readouterr()
        payload = json.loads(captured.out)
        assert np.isfinite(payload["r_squared_total"])
        assert np.isfinite(payload["r_squared_phylo"])

    def test_r2_in_text_output(self, capsys):
        """Text output includes R²_phylo and R²_total strings."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            response="brain_size",
            predictors=["body_mass"],
            method="BM",
            json=False,
        )
        svc = PhylogeneticRegression(args)
        svc.run()
        captured = capsys.readouterr()
        assert "R-squared (total):" in captured.out
        assert "R-squared (phylo):" in captured.out

    def test_r2_matches_manual_gls_reference(self, capsys):
        """R² values must match manual GLS computation (BM, no lambda estimation).

        Reference (tests/r_validation/validate_pgls_r2.R, manual GLS section):
          sigma2_ols = var(y)*n/n (ML)
          sigma2_gls_full = residuals' C_inv residuals / n
          r2_total = 1 - sigma2_gls_full / sigma2_ols
          r2_pred = existing GLS R²
          r2_phylo = r2_total - r2_pred

        Python computed values (brain_size ~ body_mass, BM):
          r_squared       = 0.9762780781
          r_squared_total = 0.9988062139
          r_squared_phylo = 0.0225281357
        """
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            response="brain_size",
            predictors=["body_mass"],
            method="BM",
            json=True,
        )
        svc = PhylogeneticRegression(args)
        svc.run()
        captured = capsys.readouterr()
        payload = json.loads(captured.out)
        assert payload["r_squared"] == pytest.approx(0.9762780781, abs=1e-6)
        assert payload["r_squared_total"] == pytest.approx(0.9988062139, abs=1e-6)
        assert payload["r_squared_phylo"] == pytest.approx(0.0225281357, abs=1e-6)
