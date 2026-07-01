import builtins
import importlib
import json
import subprocess
import sys
from io import StringIO

import pytest
import numpy as np
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch
from Bio import Phylo

from phykit.services.tree.phylo_logistic import PhyloLogistic
import phykit.services.tree.phylo_logistic as phylo_logistic_module
from phykit.helpers.trait_parsing import parse_multi_trait_file
from phykit.errors import PhykitUserError


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
BINARY_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_binary_traits.tsv")
GLM_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_glm_traits.tsv")


def test_binary_response_class_counts_uses_count_nonzero(monkeypatch):
    observed = []
    original_count_nonzero = np.count_nonzero

    def count_nonzero_spy(values, *args, **kwargs):
        observed.append(values.copy())
        return original_count_nonzero(values, *args, **kwargs)

    def fail_sum(*_args, **_kwargs):
        raise AssertionError("binary response counts should use count_nonzero")

    monkeypatch.setattr(
        phylo_logistic_module.np,
        "count_nonzero",
        count_nonzero_spy,
        raising=False,
    )
    monkeypatch.setattr(phylo_logistic_module.np, "sum", fail_sum, raising=False)

    y = np.array([0, 1, 1, 0, 1], dtype=np.int8)

    assert phylo_logistic_module._binary_response_class_counts(y) == (3, 2)
    assert len(observed) == 1
    assert observed[0].tolist() == [0, 1, 1, 0, 1]


def test_module_import_does_not_import_numpy_or_scipy():
    code = """
import sys
import phykit.services.tree.phylo_logistic as module
assert hasattr(module.np, "__getattr__")
assert module._CHO_FACTOR is None
assert module._CHO_SOLVE is None
assert module._MINIMIZE is None
assert callable(module.print_json)
assert callable(module.parse_multi_trait_file)
assert "pickle" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.trait_parsing" not in sys.modules
assert "typing" not in sys.modules
assert "numpy" not in sys.modules
assert "scipy.linalg" not in sys.modules
assert "scipy.optimize" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_module_import_does_not_import_scipy_linalg_or_optimize(monkeypatch):
    module_name = "phykit.services.tree.phylo_logistic"
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
                "phylo_logistic module import should not import SciPy linalg/optimize"
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


def test_repeated_info_matrix_caches_scipy_linalg_imports(basic_args, monkeypatch):
    previous_cho_factor = phylo_logistic_module._CHO_FACTOR
    previous_cho_solve = phylo_logistic_module._CHO_SOLVE
    phylo_logistic_module._CHO_FACTOR = None
    phylo_logistic_module._CHO_SOLVE = None
    original_import = builtins.__import__
    scipy_linalg_imports = 0

    def counting_import(name, globals=None, locals=None, fromlist=(), level=0):
        nonlocal scipy_linalg_imports
        if name == "scipy.linalg":
            scipy_linalg_imports += 1
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", counting_import)
    try:
        svc = PhyloLogistic(basic_args)
        vcv = np.eye(4)
        X = np.column_stack([np.ones(4), np.arange(4.0)])
        mu = np.full(4, 0.5)
        svc._compute_info_matrix_cholesky(X, mu, vcv)
        first_call_imports = scipy_linalg_imports
        svc._compute_info_matrix_cholesky(X, mu, vcv)
    finally:
        phylo_logistic_module._CHO_FACTOR = previous_cho_factor
        phylo_logistic_module._CHO_SOLVE = previous_cho_solve

    assert first_call_imports > 0
    assert scipy_linalg_imports == first_call_imports


def test_repeated_minimize_wrapper_caches_scipy_optimize_import(monkeypatch):
    previous_minimize = phylo_logistic_module._MINIMIZE
    phylo_logistic_module._MINIMIZE = None
    original_import = builtins.__import__
    scipy_optimize_imports = 0

    def counting_import(name, globals=None, locals=None, fromlist=(), level=0):
        nonlocal scipy_optimize_imports
        if name == "scipy.optimize":
            scipy_optimize_imports += 1
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", counting_import)
    try:
        objective = lambda x: (x[0] - 1.0) ** 2
        phylo_logistic_module.minimize(
            objective, [0.0], method="Nelder-Mead", options={"maxiter": 1}
        )
        first_call_imports = scipy_optimize_imports
        phylo_logistic_module.minimize(
            objective, [0.0], method="Nelder-Mead", options={"maxiter": 1}
        )
    finally:
        phylo_logistic_module._MINIMIZE = previous_minimize

    assert first_call_imports > 0
    assert scipy_optimize_imports == first_call_imports


def test_normal_two_tailed_p_values_match_expected_values():
    p_values = phylo_logistic_module._normal_two_tailed_p_values(
        np.array([0.0, 1.959963984540054, -2.5758293035489004])
    )
    assert p_values == pytest.approx([1.0, 0.05, 0.01])


def test_normal_two_tailed_p_values_do_not_import_scipy_stats(monkeypatch):
    original_import = __import__

    def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "scipy.stats" or name.startswith("scipy.stats."):
            raise AssertionError("phylo logistic z-test p-values should not import scipy.stats")
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr("builtins.__import__", fake_import)

    p_values = phylo_logistic_module._normal_two_tailed_p_values(
        np.array([0.0, 1.959963984540054])
    )
    assert p_values == pytest.approx([1.0, 0.05])


def test_normal_two_tailed_p_values_uses_vectorized_special_erfc(monkeypatch):
    captured = {}

    def fake_special_erfc(values):
        captured["values"] = values
        return np.full(values.shape, 0.25)

    monkeypatch.setattr(phylo_logistic_module, "special_erfc", fake_special_erfc)

    p_values = phylo_logistic_module._normal_two_tailed_p_values(
        np.array([0.0, -2.0, 4.0])
    )

    assert captured["values"].shape == (3,)
    np.testing.assert_allclose(
        captured["values"],
        np.array([0.0, 2.0, 4.0]) / np.sqrt(2.0),
    )
    np.testing.assert_allclose(p_values, np.array([0.25, 0.25, 0.25]))


def test_standard_errors_from_info_matrix_cholesky_matches_inverse():
    rng = np.random.default_rng(20260628)
    A = rng.normal(size=(8, 8))
    info_matrix = A @ A.T + np.eye(8) * 0.1

    expected = np.sqrt(np.abs(np.diag(np.linalg.inv(info_matrix))))
    observed = phylo_logistic_module._standard_errors_from_info_matrix(info_matrix)

    np.testing.assert_allclose(observed, expected)


def test_standard_errors_from_info_matrix_spd_avoids_inverse(monkeypatch):
    rng = np.random.default_rng(20260628)
    A = rng.normal(size=(5, 5))
    info_matrix = A @ A.T + np.eye(5)

    def fail_inverse(_matrix):
        raise AssertionError("SPD standard errors should use the Cholesky path")

    monkeypatch.setattr(phylo_logistic_module.np.linalg, "inv", fail_inverse)

    observed = phylo_logistic_module._standard_errors_from_info_matrix(info_matrix)

    assert np.all(np.isfinite(observed))


def test_standard_errors_from_info_matrix_keeps_inverse_fallback():
    info_matrix = np.array([[1.0, 2.0], [2.0, 1.0]])

    expected = np.sqrt(np.abs(np.diag(np.linalg.inv(info_matrix))))
    observed = phylo_logistic_module._standard_errors_from_info_matrix(info_matrix)

    np.testing.assert_allclose(observed, expected)


@pytest.fixture
def basic_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=BINARY_TRAITS_FILE,
        response="has_wings",
        predictor="body_mass",
        method="logistic_MPLE",
        json=False,
    )


@pytest.fixture
def json_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=BINARY_TRAITS_FILE,
        response="has_wings",
        predictor="body_mass",
        method="logistic_MPLE",
        json=True,
    )


@pytest.fixture
def glm_traits_args():
    """Uses the GLM traits file which has binary_trait column."""
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=GLM_TRAITS_FILE,
        response="binary_trait",
        predictor="body_mass",
        method="logistic_MPLE",
        json=False,
    )


@pytest.fixture
def multi_predictor_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=GLM_TRAITS_FILE,
        response="binary_trait",
        predictor="body_mass,count_trait",
        method="logistic_MPLE",
        json=False,
    )


class TestProcessArgs:
    def test_defaults(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            response="y",
            predictor="x1",
            method="logistic_MPLE",
        )
        svc = PhyloLogistic.__new__(PhyloLogistic)
        parsed = svc.process_args(args)
        assert parsed["tree_file_path"] == "t.tre"
        assert parsed["trait_data_path"] == "d.tsv"
        assert parsed["response"] == "y"
        assert parsed["predictors"] == ["x1"]
        assert parsed["method"] == "logistic_MPLE"
        assert parsed["json_output"] is False

    def test_comma_separated_predictors(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            response="y",
            predictor="x1,x2,x3",
            method="logistic_MPLE",
        )
        svc = PhyloLogistic.__new__(PhyloLogistic)
        parsed = svc.process_args(args)
        assert parsed["predictors"] == ["x1", "x2", "x3"]


class TestTransformedVCV:
    def test_transformed_vcv_alpha_zero(self, basic_args):
        """When alpha ~ 0, the transformed VCV should be the standard BM VCV."""
        from phykit.services.tree.vcv_utils import build_vcv_matrix

        svc = PhyloLogistic(basic_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(BINARY_TRAITS_FILE, tree_tips)
        ordered_names = sorted(traits.keys())

        vcv_bm = build_vcv_matrix(tree, ordered_names)
        vcv_ou, diag_corr = svc._build_logistic_vcv(
            tree, 1e-12, ordered_names, build_vcv_matrix
        )

        # diag_corr should be all ones for alpha ~ 0
        np.testing.assert_allclose(diag_corr, np.ones(len(ordered_names)), atol=1e-6)
        # VCV should be the standard BM VCV
        np.testing.assert_allclose(vcv_ou, vcv_bm, atol=1e-4)

    def test_transformed_vcv_positive_alpha(self, basic_args):
        """Positive alpha should give a valid (PD) VCV matrix."""
        from phykit.services.tree.vcv_utils import build_vcv_matrix

        svc = PhyloLogistic(basic_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(BINARY_TRAITS_FILE, tree_tips)
        ordered_names = sorted(traits.keys())

        vcv, diag_corr = svc._build_logistic_vcv(
            tree, 0.05, ordered_names, build_vcv_matrix
        )

        # VCV should be symmetric
        np.testing.assert_allclose(vcv, vcv.T, atol=1e-10)

        # VCV should be positive definite
        eigenvalues = np.linalg.eigvalsh(vcv)
        assert np.all(eigenvalues > -1e-10), f"VCV has negative eigenvalue: {min(eigenvalues)}"

        # Diagonal correction should be positive and < 1
        assert np.all(diag_corr > 0)
        assert np.all(diag_corr <= 1.0 + 1e-10)

    def test_transformed_vcv_matches_copied_tree_transform(self, basic_args):
        """Direct transformed VCV should match the previous copy/mutate path."""
        import pickle
        from phykit.services.tree.vcv_utils import build_vcv_matrix

        svc = PhyloLogistic(basic_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        ordered_names = ["A", "B", "C", "D"]
        alpha = 0.05

        expected_tree = pickle.loads(
            pickle.dumps(tree, protocol=pickle.HIGHEST_PROTOCOL)
        )
        for clade in expected_tree.find_clades():
            if clade.branch_length is not None and clade.branch_length > 0:
                bl = clade.branch_length
                clade.branch_length = (1 - np.exp(-2 * alpha * bl)) / (2 * alpha)
        expected = build_vcv_matrix(expected_tree, ordered_names)
        diag_corr = np.exp(
            -2 * alpha * svc._root_tip_distances(tree, ordered_names)
        )
        np.fill_diagonal(expected, np.diag(expected) + diag_corr)

        original_lengths = [
            clade.branch_length for clade in tree.find_clades(order="preorder")
        ]
        observed, _ = svc._build_logistic_vcv(
            tree, alpha, ordered_names, build_vcv_matrix
        )
        after_lengths = [
            clade.branch_length for clade in tree.find_clades(order="preorder")
        ]

        np.testing.assert_allclose(observed, expected)
        assert after_lengths == original_lengths

    def test_root_tip_distances_fast_path_does_not_call_distance(
        self, basic_args, monkeypatch
    ):
        svc = PhyloLogistic(basic_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,C:2):0.5;"), "newick")

        def fail_depths(*args, **kwargs):
            raise AssertionError("depths fallback should not be called")

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("terminal fallback should not be called")

        def fail_distance(*args, **kwargs):
            raise AssertionError("distance fallback should not be called")

        monkeypatch.setattr(tree, "depths", fail_depths)
        monkeypatch.setattr(tree, "get_terminals", fail_get_terminals)
        monkeypatch.setattr(tree, "distance", fail_distance)
        distances = svc._root_tip_distances(tree, ["A", "B", "C"])

        np.testing.assert_allclose(distances, np.array([2.0, 2.0, 2.0]))

    def test_fit_reuses_root_tip_distances_for_vcv_builds(
        self, basic_args, monkeypatch
    ):
        from phykit.services.tree.vcv_utils import build_vcv_matrix

        svc = PhyloLogistic(basic_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        ordered_names = ["A", "B", "C", "D"]
        X = np.column_stack([np.ones(4), np.array([0.1, 0.2, 0.8, 0.9])])
        y = np.array([0.0, 0.0, 1.0, 1.0])
        call_count = 0
        original = PhyloLogistic._root_tip_distances

        def counting_root_tip_distances(counted_tree, counted_names):
            nonlocal call_count
            call_count += 1
            return original(counted_tree, counted_names)

        def fake_minimize(_objective, x0, *args, **kwargs):
            return Namespace(x=x0, status=0)

        monkeypatch.setattr(
            PhyloLogistic,
            "_root_tip_distances",
            staticmethod(counting_root_tip_distances),
        )
        monkeypatch.setattr(phylo_logistic_module, "minimize", fake_minimize)

        svc._fit(tree, y, X, ordered_names, build_vcv_matrix)

        assert call_count == 1

    def test_logistic_vcv_diag_correction_fast_path_does_not_call_distance(
        self, basic_args, monkeypatch
    ):
        from phykit.services.tree.vcv_utils import build_vcv_matrix

        svc = PhyloLogistic(basic_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,C:2):0.5;"), "newick")

        def fail_distance(self, *args, **kwargs):
            raise AssertionError("distance fallback should not be called")

        monkeypatch.setattr(type(tree), "distance", fail_distance)
        _, diag_corr = svc._build_logistic_vcv(
            tree, 0.05, ["A", "B", "C"], build_vcv_matrix
        )

        expected = np.exp(-2 * 0.05 * np.array([2.0, 2.0, 2.0]))
        np.testing.assert_allclose(diag_corr, expected)

    def test_info_matrix_cholesky_matches_inverse(self, basic_args):
        rng = np.random.default_rng(20260623)
        svc = PhyloLogistic(basic_args)
        n = 16
        A = rng.normal(size=(n, n))
        vcv = A @ A.T + np.eye(n)
        X = np.column_stack([np.ones(n), rng.normal(size=(n, 2))])
        eta = rng.normal(size=n) * 0.25
        mu = 1.0 / (1.0 + np.exp(-eta))

        fast = svc._compute_info_matrix_cholesky(X, mu, vcv)
        inverse = svc._compute_info_matrix_inverse(X, mu, vcv)

        assert fast == pytest.approx(inverse)

    def test_info_matrix_inverse_matches_diagonal_weight_reference(self, basic_args):
        rng = np.random.default_rng(20260623)
        svc = PhyloLogistic(basic_args)
        n = 24
        A = rng.normal(size=(n, n))
        vcv = A @ A.T + np.eye(n)
        X = np.column_stack([np.ones(n), rng.normal(size=(n, 4))])
        eta = rng.normal(size=n) * 0.25
        mu = 1.0 / (1.0 + np.exp(-eta))

        W_sqrt = np.diag(np.sqrt(mu * (1 - mu)))
        expected = X.T @ W_sqrt @ np.linalg.inv(vcv) @ W_sqrt @ X
        observed = svc._compute_info_matrix_inverse(X, mu, vcv)

        np.testing.assert_allclose(observed, expected)


class TestFit:
    def test_logistic_starting_values_match_diagonal_weight_reference(self):
        svc = PhyloLogistic.__new__(PhyloLogistic)
        rng = np.random.default_rng(20260623)
        n = 40
        X = np.column_stack([np.ones(n), rng.normal(size=(n, 3))])
        y = rng.binomial(1, 0.45, size=n).astype(float)

        def reference():
            beta = np.zeros(X.shape[1])
            for _ in range(50):
                eta = np.clip(X @ beta, -10, 10)
                mu = 1.0 / (1.0 + np.exp(-eta))
                mu = np.clip(mu, 1e-6, 1 - 1e-6)
                weights = np.diag(mu * (1 - mu))
                z = eta + (y - mu) / (mu * (1 - mu))
                beta_new = np.linalg.solve(X.T @ weights @ X, X.T @ weights @ z)
                if np.max(np.abs(beta_new - beta)) < 1e-8:
                    return beta_new
                beta = beta_new
            return beta

        observed = svc._logistic_starting_values(y, X, btol=10)
        np.testing.assert_allclose(observed, reference())

    def test_fit_returns_coefficients(self, glm_traits_args):
        """Fitting should return finite beta, SE, z, p for each coefficient."""
        from phykit.services.tree.vcv_utils import build_vcv_matrix

        svc = PhyloLogistic(glm_traits_args)
        tree = svc.read_tree_file()
        svc.validate_tree(tree, min_tips=3, require_branch_lengths=True)
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(GLM_TRAITS_FILE, tree_tips)
        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        resp_idx = trait_names.index("binary_trait")
        pred_idx = trait_names.index("body_mass")
        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X = np.column_stack([
            np.ones(n),
            np.array([traits[name][pred_idx] for name in ordered_names]),
        ])

        result = svc._fit(tree, y, X, ordered_names, build_vcv_matrix)

        assert "coefficients" in result
        intercept = result["coefficients"]["(Intercept)"]
        assert np.isfinite(intercept["estimate"])
        assert np.isfinite(intercept["std_error"])
        assert intercept["std_error"] > 0
        assert np.isfinite(intercept["z_value"])
        assert np.isfinite(intercept["p_value"])
        assert 0 <= intercept["p_value"] <= 1

    def test_alpha_positive(self, glm_traits_args):
        """Estimated alpha should be positive."""
        from phykit.services.tree.vcv_utils import build_vcv_matrix

        svc = PhyloLogistic(glm_traits_args)
        tree = svc.read_tree_file()
        svc.validate_tree(tree, min_tips=3, require_branch_lengths=True)
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(GLM_TRAITS_FILE, tree_tips)
        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        resp_idx = trait_names.index("binary_trait")
        pred_idx = trait_names.index("body_mass")
        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X = np.column_stack([
            np.ones(n),
            np.array([traits[name][pred_idx] for name in ordered_names]),
        ])

        result = svc._fit(tree, y, X, ordered_names, build_vcv_matrix)
        assert result["alpha"] > 0

    def test_convergence(self, glm_traits_args):
        """Optimizer should converge (code 0 or 2 for L-BFGS-B)."""
        from phykit.services.tree.vcv_utils import build_vcv_matrix

        svc = PhyloLogistic(glm_traits_args)
        tree = svc.read_tree_file()
        svc.validate_tree(tree, min_tips=3, require_branch_lengths=True)
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(GLM_TRAITS_FILE, tree_tips)
        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        resp_idx = trait_names.index("binary_trait")
        pred_idx = trait_names.index("body_mass")
        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X = np.column_stack([
            np.ones(n),
            np.array([traits[name][pred_idx] for name in ordered_names]),
        ])

        result = svc._fit(tree, y, X, ordered_names, build_vcv_matrix)
        # 0 = converged, 1 = max iterations, 2 = abnormal termination in line search
        # For small datasets, code 2 is acceptable (numerical precision at boundary)
        assert result["convergence"] in (0, 2)


class TestValidation:
    def test_missing_column(self, basic_args):
        basic_args.response = "nonexistent_column"
        svc = PhyloLogistic(basic_args)
        with pytest.raises(PhykitUserError) as exc_info:
            svc.run()
        assert any("nonexistent_column" in m for m in exc_info.value.messages)

    def test_response_in_predictors(self):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=BINARY_TRAITS_FILE,
            response="has_wings",
            predictor="has_wings",
            method="logistic_MPLE",
            json=False,
        )
        svc = PhyloLogistic(args)
        with pytest.raises(PhykitUserError) as exc_info:
            svc.run()
        assert any("must not also be a predictor" in m for m in exc_info.value.messages)


class TestRun:
    def test_run_uses_unmodified_tree_read(self, mocker):
        args = Namespace(
            tree="/some/path/to/file.tre",
            trait_data="/some/path/to/traits.tsv",
            response="y",
            predictor="x",
            method="logistic_MPLE",
            json=False,
        )
        svc = PhyloLogistic(args)
        tree = object()
        mocked_read = mocker.patch.object(
            PhyloLogistic,
            "read_tree_file_unmodified",
            return_value=tree,
        )
        mocker.patch.object(
            PhyloLogistic,
            "read_tree_file",
            side_effect=AssertionError("run should not copy cached trees"),
        )
        mocked_validate = mocker.patch.object(
            PhyloLogistic, "validate_tree"
        )
        mocker.patch.object(
            PhyloLogistic,
            "get_tip_names_from_tree",
            return_value=["a", "b", "c"],
        )
        mocker.patch(
            "phykit.services.tree.phylo_logistic.parse_multi_trait_file",
            return_value=(
                ["y", "x"],
                {
                    "a": [0.0, 1.0],
                    "b": [1.0, 2.0],
                    "c": [0.0, 3.0],
                },
            ),
        )
        fit_result = {
            "coefficients": {},
            "log_likelihood": -1.0,
            "alpha": 1.0,
            "convergence": 0,
        }
        mocked_fit = mocker.patch.object(
            PhyloLogistic,
            "_fit",
            return_value=fit_result,
        )
        mocker.patch.object(PhyloLogistic, "_print_text_output")

        svc.run()

        mocked_read.assert_called_once_with()
        mocked_validate.assert_called_once_with(
            tree,
            min_tips=3,
            require_branch_lengths=True,
            context="phylogenetic logistic regression",
        )
        assert mocked_fit.call_args.args[0] is tree

    def test_run_resolves_trait_columns_with_first_index_map(self, mocker):
        args = Namespace(
            tree="/some/path/to/file.tre",
            trait_data="/some/path/to/traits.tsv",
            response="y",
            predictor="x",
            method="logistic_MPLE",
            json=False,
        )
        svc = PhyloLogistic(args)
        tree = object()
        mocker.patch.object(svc, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(svc, "validate_tree")
        mocker.patch.object(svc, "get_tip_names_from_tree", return_value=["a", "b", "c"])
        mocker.patch(
            "phykit.services.tree.phylo_logistic.parse_multi_trait_file",
            return_value=(
                ["y", "x", "y", "x"],
                {
                    "a": [0.0, 1.0, 9.0, 9.0],
                    "b": [1.0, 2.0, 9.0, 9.0],
                    "c": [0.0, 3.0, 9.0, 9.0],
                },
            ),
        )

        def arrays(traits, ordered_names, response_index, predictor_indices):
            assert response_index == 0
            assert predictor_indices == [1]
            return np.array([0.0, 1.0, 0.0]), np.ones((3, 2))

        mocker.patch(
            "phykit.services.tree.phylo_logistic.response_predictor_arrays",
            side_effect=arrays,
        )
        mocker.patch.object(
            PhyloLogistic,
            "_fit",
            return_value={
                "coefficients": {},
                "log_likelihood": -1.0,
                "alpha": 1.0,
                "convergence": 0,
            },
        )
        mocker.patch.object(PhyloLogistic, "_print_text_output")

        svc.run()

    @patch("builtins.print")
    def test_creates_results(self, mocked_print, glm_traits_args):
        """Text output should contain key sections."""
        svc = PhyloLogistic(glm_traits_args)
        svc.run()
        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Phylogenetic Logistic Regression" in all_output
        assert "(Intercept)" in all_output
        assert "body_mass" in all_output
        assert "Alpha:" in all_output
        assert "Log-likelihood:" in all_output
        assert "Penalized log-likelihood:" in all_output
        assert "AIC:" in all_output
        assert "Signif. codes:" in all_output

    @patch("builtins.print")
    def test_print_text_output_batches_coefficient_rows(self, mocked_print):
        svc = PhyloLogistic.__new__(PhyloLogistic)
        svc.response = "binary_trait"
        svc.predictors = ["body_mass", "count_trait"]
        result = {
            "formula": "binary_trait ~ body_mass + count_trait",
            "method": "logistic_MPLE",
            "n_taxa": 25,
            "alpha": 0.125,
            "coefficients": {
                "(Intercept)": {
                    "estimate": 1.0,
                    "std_error": 0.25,
                    "z_value": 4.0,
                    "p_value": 0.0005,
                },
                "body_mass": {
                    "estimate": -0.5,
                    "std_error": 0.2,
                    "z_value": -2.5,
                    "p_value": 0.012,
                },
            },
            "log_likelihood": -12.5,
            "penalized_log_likelihood": -13.5,
            "aic": 31.0,
        }

        svc._print_text_output(result)

        mocked_print.assert_called_once_with(
            "Phylogenetic Logistic Regression (Ives & Garland 2010)\n"
            "======================================================\n"
            "Response: binary_trait\n"
            "Predictor(s): body_mass, count_trait\n"
            "Formula: binary_trait ~ body_mass + count_trait\n"
            "Method: logistic_MPLE\n"
            "Taxa: 25\n"
            "Alpha: 0.125000\n"
            "\n"
            "Coefficients:\n"
            "                        Estimate   Std.Error     z-value    Pr(>|z|)\n"
            "(Intercept)               1.0000      0.2500      4.0000    0.000500    ***\n"
            "body_mass                -0.5000      0.2000     -2.5000    0.012000    *\n"
            "---\n"
            "Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1\n"
            "\n"
            "Log-likelihood: -12.5000\n"
            "Penalized log-likelihood: -13.5000\n"
            "AIC: 31.0000"
        )

    def test_json_output(self, capsys):
        """JSON should have all expected fields."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=GLM_TRAITS_FILE,
            response="binary_trait",
            predictor="body_mass",
            method="logistic_MPLE",
            json=True,
        )
        svc = PhyloLogistic(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert data["method"] == "logistic_MPLE"
        assert data["response"] == "binary_trait"
        assert data["predictors"] == ["body_mass"]
        assert "alpha" in data
        assert data["alpha"] > 0
        assert "(Intercept)" in data["coefficients"]
        assert "body_mass" in data["coefficients"]
        assert "log_likelihood" in data
        assert "penalized_log_likelihood" in data
        assert "aic" in data
        assert "convergence" in data
        assert "n_taxa" in data

    def test_format_result_maps_coefficients_with_zip(self):
        svc = PhyloLogistic.__new__(PhyloLogistic)
        svc.response = "binary_trait"
        svc.predictors = ["body_mass"]

        result = svc._format_result(
            method="logistic_MPLE",
            coef_names=["(Intercept)", "body_mass"],
            beta_hat=np.array([1.25, -0.5]),
            se=np.array([0.25, 0.2]),
            z_stats=np.array([5.0, -2.5]),
            p_values=np.array([0.0001, 0.012]),
            ll=-12.5,
            pen_ll=-13.5,
            aic=31.0,
            formula="binary_trait ~ body_mass",
            n=25,
            k=1,
            ordered_names=["a", "b"],
            alpha=0.125,
            convergence=0,
        )

        assert result["coefficients"] == {
            "(Intercept)": {
                "estimate": 1.25,
                "std_error": 0.25,
                "z_value": 5.0,
                "p_value": 0.0001,
            },
            "body_mass": {
                "estimate": -0.5,
                "std_error": 0.2,
                "z_value": -2.5,
                "p_value": 0.012,
            },
        }

    @patch("builtins.print")
    def test_multiple_predictors(self, mocked_print, multi_predictor_args):
        """Comma-separated predictors should work."""
        svc = PhyloLogistic(multi_predictor_args)
        svc.run()
        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "body_mass" in all_output
        assert "count_trait" in all_output
        assert "binary_trait ~ body_mass + count_trait" in all_output
