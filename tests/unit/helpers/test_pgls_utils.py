import builtins
import importlib
import sys
from io import StringIO

import numpy as np
import pytest
from Bio import Phylo

from phykit.helpers.pgls_utils import (
    _max_lambda_fallback,
    _pgls_log_likelihood_inverse,
    estimate_lambda,
    fit_gls,
    max_lambda,
    pgls_log_likelihood,
)
from phykit.errors import PhykitUserError


def test_module_import_does_not_import_scipy_linalg_or_optimize(monkeypatch):
    module_name = "phykit.helpers.pgls_utils"
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
                "pgls_utils module import should not import SciPy linalg/optimize"
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


def test_module_import_does_not_import_numpy_or_scipy():
    code = """
import sys
import phykit.helpers.pgls_utils as module

assert hasattr(module.np, "__getattr__")
assert "typing" not in sys.modules
assert "numpy" not in sys.modules
assert "scipy" not in sys.modules
assert "scipy.linalg" not in sys.modules
assert "scipy.optimize" not in sys.modules
"""
    import subprocess

    subprocess.run([sys.executable, "-c", code], check=True)


def test_max_lambda_fast_path_matches_fallback_without_distance(monkeypatch):
    tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
    expected = _max_lambda_fallback(tree)

    def fail_distance(*args, **kwargs):
        raise AssertionError("distance should not be called")

    monkeypatch.setattr(tree, "distance", fail_distance)

    assert max_lambda(tree) == pytest.approx(expected)


def test_max_lambda_fast_path_detects_non_ultrametric_without_distance(monkeypatch):
    tree = Phylo.read(StringIO("((A:1,B:2):1,(C:1,D:1):1);"), "newick")

    def fail_distance(*args, **kwargs):
        raise AssertionError("distance should not be called")

    monkeypatch.setattr(tree, "distance", fail_distance)

    assert max_lambda(tree) == 1.0


def test_max_lambda_uses_direct_tree_traversal(monkeypatch):
    tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")

    def fail_traversal(*_args, **_kwargs):
        raise AssertionError("generic tree traversal should not be used")

    monkeypatch.setattr(tree, "get_terminals", fail_traversal)
    monkeypatch.setattr(tree, "depths", fail_traversal)
    monkeypatch.setattr(tree, "find_clades", fail_traversal)
    monkeypatch.setattr(tree, "distance", fail_traversal)

    assert max_lambda(tree) == pytest.approx(2.0)


def test_max_lambda_direct_path_handles_mixed_child_counts(monkeypatch):
    tree = Phylo.read(
        StringIO("(A:2,(B:1,C:1):1,(D:1,E:1,F:1):1);"),
        "newick",
    )

    def fail_traversal(*_args, **_kwargs):
        raise AssertionError("generic tree traversal should not be used")

    monkeypatch.setattr(tree, "get_terminals", fail_traversal)
    monkeypatch.setattr(tree, "depths", fail_traversal)
    monkeypatch.setattr(tree, "find_clades", fail_traversal)
    monkeypatch.setattr(tree, "distance", fail_traversal)

    assert max_lambda(tree) == pytest.approx(2.0)


def test_pgls_log_likelihood_cholesky_matches_inverse_formula():
    rng = np.random.default_rng(20260622)
    n = 40
    A = rng.normal(size=(n, n))
    C = A @ A.T + np.eye(n) * 0.5
    X = np.column_stack([np.ones(n), rng.normal(size=(n, 2))])
    y = X @ np.array([0.5, 1.2, -0.7]) + rng.normal(scale=0.1, size=n)

    assert pgls_log_likelihood(y, X, C) == pytest.approx(
        _pgls_log_likelihood_inverse(y, X, C),
        abs=1e-8,
    )


def test_pgls_log_likelihood_uses_single_cholesky_solve(monkeypatch):
    import phykit.helpers.pgls_utils as pgls_utils

    rng = np.random.default_rng(20260628)
    n = 40
    A = rng.normal(size=(n, n))
    C = A @ A.T + np.eye(n) * 0.5
    X = np.column_stack([np.ones(n), rng.normal(size=(n, 2))])
    y = X @ np.array([0.5, 1.2, -0.7]) + rng.normal(scale=0.1, size=n)
    original_cho_solve = pgls_utils.cho_solve
    solve_calls = 0

    def counting_cho_solve(*args, **kwargs):
        nonlocal solve_calls
        solve_calls += 1
        return original_cho_solve(*args, **kwargs)

    monkeypatch.setattr(pgls_utils, "cho_solve", counting_cho_solve)

    observed = pgls_log_likelihood(y, X, C)
    expected = _pgls_log_likelihood_inverse(y, X, C)

    assert solve_calls == 1
    assert observed == pytest.approx(expected, abs=1e-8)


def test_pgls_log_likelihood_singular_matrix_returns_floor():
    y = np.array([1.0, 2.0, 3.0])
    X = np.column_stack([np.ones(3), np.array([0.0, 1.0, 2.0])])
    C = np.ones((3, 3))

    assert pgls_log_likelihood(y, X, C) == -1e20


def test_repeated_pgls_log_likelihood_caches_scipy_linalg_imports(monkeypatch):
    import phykit.helpers.pgls_utils as pgls_utils

    previous_cho_factor = pgls_utils._CHO_FACTOR
    previous_cho_solve = pgls_utils._CHO_SOLVE
    pgls_utils._CHO_FACTOR = None
    pgls_utils._CHO_SOLVE = None
    original_import = builtins.__import__
    scipy_linalg_imports = 0

    def counting_import(name, globals=None, locals=None, fromlist=(), level=0):
        nonlocal scipy_linalg_imports
        if name == "scipy.linalg":
            scipy_linalg_imports += 1
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", counting_import)
    try:
        C = np.eye(4)
        X = np.column_stack([np.ones(4), np.arange(4.0)])
        y = np.array([1.0, 1.5, 2.0, 2.5])
        pgls_log_likelihood(y, X, C)
        pgls_log_likelihood(y, X, C)
    finally:
        pgls_utils._CHO_FACTOR = previous_cho_factor
        pgls_utils._CHO_SOLVE = previous_cho_solve

    assert scipy_linalg_imports == 2


def test_repeated_minimize_scalar_wrapper_caches_scipy_optimize_import(monkeypatch):
    import phykit.helpers.pgls_utils as pgls_utils

    previous_minimize_scalar = pgls_utils._MINIMIZE_SCALAR
    pgls_utils._MINIMIZE_SCALAR = None
    original_import = builtins.__import__
    scipy_optimize_imports = 0

    def counting_import(name, globals=None, locals=None, fromlist=(), level=0):
        nonlocal scipy_optimize_imports
        if name == "scipy.optimize":
            scipy_optimize_imports += 1
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", counting_import)
    try:
        objective = lambda x: (x - 1.0) ** 2
        pgls_utils.minimize_scalar(
            objective, bounds=(0.0, 2.0), method="bounded"
        )
        first_call_imports = scipy_optimize_imports
        pgls_utils.minimize_scalar(
            objective, bounds=(0.0, 2.0), method="bounded"
        )
    finally:
        pgls_utils._MINIMIZE_SCALAR = previous_minimize_scalar

    assert first_call_imports > 0
    assert scipy_optimize_imports == first_call_imports


def test_estimate_lambda_restores_diagonal_directly(monkeypatch):
    import phykit.helpers.pgls_utils as pgls_utils

    vcv = np.array(
        [
            [2.0, 0.4, 0.2],
            [0.4, 3.0, 0.5],
            [0.2, 0.5, 4.0],
        ]
    )
    y = np.array([1.0, 1.5, 2.0])
    X = np.column_stack([np.ones(3), np.arange(3.0)])
    original_diag = vcv.diagonal().copy()
    likelihood_calls = 0

    def fail_diag(*_args, **_kwargs):
        raise AssertionError("estimate_lambda should use ndarray diagonal access")

    def fail_fill_diagonal(*_args, **_kwargs):
        raise AssertionError("estimate_lambda should restore diagonal directly")

    def fake_likelihood(_y, _X, C_lam):
        nonlocal likelihood_calls
        likelihood_calls += 1
        np.testing.assert_allclose(C_lam.diagonal(), original_diag)
        return -abs(C_lam[0, 1] - 0.2)

    def fake_minimize_scalar(fn, bounds, method):
        assert method == "bounded"
        x = sum(bounds) / 2.0
        return type("Result", (), {"x": x, "fun": fn(x)})()

    monkeypatch.setattr(pgls_utils.np, "diag", fail_diag)
    monkeypatch.setattr(pgls_utils.np, "fill_diagonal", fail_fill_diagonal)
    monkeypatch.setattr(pgls_utils, "pgls_log_likelihood", fake_likelihood)
    monkeypatch.setattr(pgls_utils, "minimize_scalar", fake_minimize_scalar)

    lambda_val, ll = estimate_lambda(y, X, vcv, max_lam=1.0)

    assert 0.0 <= lambda_val <= 1.0
    assert np.isfinite(ll)
    assert likelihood_calls == 11


def test_fit_gls_matches_scalar_reference():
    rng = np.random.default_rng(20260625)
    n = 30
    A = rng.normal(size=(n, n))
    C_inv = np.linalg.inv(A @ A.T + np.eye(n))
    X = np.column_stack([np.ones(n), rng.normal(size=(n, 3))])
    beta = np.array([0.2, 1.0, -0.4, 0.7])
    y = X @ beta + rng.normal(scale=0.1, size=n)

    beta_hat, residuals, sigma2, var_beta = fit_gls(y, X, C_inv)

    XtCiX = X.T @ C_inv @ X
    XtCiX_inv = np.linalg.inv(XtCiX)
    expected_beta = XtCiX_inv @ (X.T @ (C_inv @ y))
    expected_residuals = y - X @ expected_beta
    df_resid = n - X.shape[1]
    expected_sigma2 = (
        float(expected_residuals @ C_inv @ expected_residuals) / df_resid
    )
    expected_var_beta = expected_sigma2 * XtCiX_inv

    assert beta_hat == pytest.approx(expected_beta)
    assert residuals == pytest.approx(expected_residuals)
    assert sigma2 == pytest.approx(expected_sigma2)
    assert var_beta == pytest.approx(expected_var_beta)


def test_fit_gls_preallocates_combined_rhs(monkeypatch):
    import phykit.helpers.pgls_utils as pgls_utils

    rng = np.random.default_rng(20260629)
    n = 18
    k = 3
    A = rng.normal(size=(n, n))
    C_inv = np.linalg.inv(A @ A.T + np.eye(n))
    X = np.empty((n, k + 1), dtype=float)
    X[:, 0] = 1.0
    X[:, 1:] = rng.normal(size=(n, k))
    beta = np.array([0.2, -0.5, 1.1, 0.7])
    y = X @ beta + rng.normal(scale=0.1, size=n)

    def fail_column_stack(*_args, **_kwargs):
        raise AssertionError("fit_gls should preallocate the combined RHS")

    monkeypatch.setattr(pgls_utils.np, "column_stack", fail_column_stack)

    beta_hat, residuals, sigma2, var_beta = fit_gls(y, X, C_inv)

    XtCiX = X.T @ C_inv @ X
    XtCiX_inv = np.linalg.inv(XtCiX)
    expected_beta = XtCiX_inv @ (X.T @ (C_inv @ y))
    expected_residuals = y - X @ expected_beta
    expected_sigma2 = (
        float(expected_residuals @ C_inv @ expected_residuals)
        / (n - X.shape[1])
    )

    assert beta_hat == pytest.approx(expected_beta)
    assert residuals == pytest.approx(expected_residuals)
    assert sigma2 == pytest.approx(expected_sigma2)
    assert var_beta == pytest.approx(expected_sigma2 * XtCiX_inv)


def test_fit_gls_singular_design_matrix_raises_user_error():
    C_inv = np.eye(4)
    X = np.column_stack([np.ones(4), np.ones(4)])
    y = np.array([1.0, 2.0, 3.0, 4.0])

    with pytest.raises(PhykitUserError) as exc_info:
        fit_gls(y, X, C_inv)

    assert "Singular design matrix: cannot estimate coefficients." in (
        exc_info.value.messages
    )
