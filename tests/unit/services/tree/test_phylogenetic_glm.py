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

from phykit.services.tree.phylogenetic_glm import PhylogeneticGLM
import phykit.services.tree.phylogenetic_glm as phylogenetic_glm_module
from phykit.helpers.trait_parsing import parse_multi_trait_file
from phykit.errors import PhykitUserError


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
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
        phylogenetic_glm_module.np,
        "count_nonzero",
        count_nonzero_spy,
        raising=False,
    )
    monkeypatch.setattr(phylogenetic_glm_module.np, "sum", fail_sum, raising=False)

    y = np.array([0, 1, 1, 0, 1], dtype=np.int8)

    assert phylogenetic_glm_module._binary_response_class_counts(y) == (3, 2)
    assert len(observed) == 1
    assert observed[0].tolist() == [0, 1, 1, 0, 1]


def test_poisson_overdispersion_uses_direct_array_sum(monkeypatch):
    y = np.array([1.0, 2.0, 4.0, 3.0, 7.0])
    mu = np.array([1.2, 1.8, 3.5, 4.2, 6.6])
    pearson_resid = (y - mu) / np.sqrt(mu)
    expected = float(np.sum(pearson_resid**2)) / 3

    def fail_sum(*_args, **_kwargs):
        raise AssertionError("Poisson overdispersion should use ndarray.sum")

    monkeypatch.setattr(phylogenetic_glm_module.np, "sum", fail_sum, raising=False)

    observed = phylogenetic_glm_module._poisson_overdispersion(y, mu, 3)

    assert observed == pytest.approx(expected)


def test_poisson_log_likelihood_uses_direct_array_sum(monkeypatch):
    from scipy.special import gammaln

    y = np.arange(260, dtype=float) % 7
    mu = np.linspace(0.5, 8.5, y.size)
    terms = y * np.log(np.clip(mu, 1e-300, None)) - mu - gammaln(y + 1)
    expected = float(terms.sum())

    def fail_sum(*_args, **_kwargs):
        raise AssertionError("Poisson likelihood should use ndarray.sum")

    monkeypatch.setattr(phylogenetic_glm_module.np, "sum", fail_sum, raising=False)

    observed = phylogenetic_glm_module._poisson_log_likelihood(y, mu)

    assert observed == pytest.approx(expected)


def test_poisson_log_likelihood_preserves_large_vector_sum_path(monkeypatch):
    y = np.arange(20001, dtype=float) % 5
    mu = np.linspace(0.5, 12.0, y.size)
    original_sum = phylogenetic_glm_module.np.sum
    calls = []

    def sum_spy(values, *args, **kwargs):
        calls.append(np.asarray(values).shape)
        return original_sum(values, *args, **kwargs)

    monkeypatch.setattr(phylogenetic_glm_module.np, "sum", sum_spy, raising=False)

    observed = phylogenetic_glm_module._poisson_log_likelihood(y, mu)

    assert np.isfinite(observed)
    assert calls == [(20001,)]


def test_ordered_distance_array_preserves_large_list_path(monkeypatch):
    count = phylogenetic_glm_module._ROOT_DISTANCE_FROMITER_MAX_TIPS + 1
    ordered_names = [f"taxon_{idx}" for idx in range(count)]
    distance_by_name = {name: float(idx) for idx, name in enumerate(ordered_names)}

    def fail_fromiter(*_args, **_kwargs):
        raise AssertionError("large distance vectors should preserve np.array path")

    monkeypatch.setattr(phylogenetic_glm_module.np, "fromiter", fail_fromiter)

    observed = phylogenetic_glm_module._ordered_distance_array(
        distance_by_name,
        ordered_names,
    )

    assert observed.shape == (count,)
    assert observed[-1] == pytest.approx(float(count - 1))


def test_mean_1d_uses_array_method_for_small_vectors(monkeypatch):
    values = np.array([1.0, 2.0, 6.0])

    def fail_mean(*_args, **_kwargs):
        raise AssertionError("small GLM vectors should use ndarray.mean")

    monkeypatch.setattr(phylogenetic_glm_module.np, "mean", fail_mean)

    assert phylogenetic_glm_module._mean_1d(values) == pytest.approx(3.0)


def test_mean_1d_preserves_large_numpy_path(monkeypatch):
    values = np.ones(phylogenetic_glm_module._DIRECT_MEAN_MAX_SIZE + 1)
    original_mean = phylogenetic_glm_module.np.mean
    calls = []

    def mean_spy(observed, *args, **kwargs):
        calls.append(observed.shape)
        return original_mean(observed, *args, **kwargs)

    monkeypatch.setattr(phylogenetic_glm_module.np, "mean", mean_spy)

    assert phylogenetic_glm_module._mean_1d(values) == pytest.approx(1.0)
    assert calls == [(phylogenetic_glm_module._DIRECT_MEAN_MAX_SIZE + 1,)]


def test_module_import_does_not_import_scipy_optimize(monkeypatch):
    module_name = "phykit.services.tree.phylogenetic_glm"
    previous = sys.modules.pop(module_name, None)
    original_import = builtins.__import__

    def guarded_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "scipy.optimize" or name.startswith("scipy.optimize."):
            raise AssertionError(
                "phylogenetic_glm module import should not import SciPy optimize"
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


def test_module_import_does_not_import_numpy_or_scipy_optimize():
    code = """
import sys
import phykit.services.tree.phylogenetic_glm as module
assert callable(module.print_json)
assert callable(module.parse_multi_trait_file)
assert hasattr(module.np, "__getattr__")
assert "typing" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.trait_parsing" not in sys.modules
assert "numpy" not in sys.modules
assert "scipy.optimize" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_lazy_numpy_caches_resolved_attributes():
    lazy_np = phylogenetic_glm_module._LazyNumpy()

    array_attr = lazy_np.array

    assert lazy_np.__dict__["array"] is array_attr
    assert lazy_np.array is array_attr
    assert lazy_np._module is not None


def test_cholesky_wrappers_cache_scipy_callables(monkeypatch):
    previous_cho_factor = phylogenetic_glm_module._CHO_FACTOR
    previous_cho_solve = phylogenetic_glm_module._CHO_SOLVE
    phylogenetic_glm_module._CHO_FACTOR = None
    phylogenetic_glm_module._CHO_SOLVE = None
    matrix = np.array([[2.0, 0.2], [0.2, 1.5]])
    rhs = np.eye(2)

    try:
        factor = phylogenetic_glm_module.cho_factor(
            matrix,
            lower=True,
            check_finite=False,
        )
        phylogenetic_glm_module.cho_solve(
            factor,
            rhs,
            check_finite=False,
        )
        original_import = builtins.__import__

        def fail_linalg_import(name, *args, **kwargs):
            if name == "scipy.linalg" or name.startswith("scipy.linalg."):
                raise AssertionError("cached wrappers should not re-import scipy.linalg")
            return original_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", fail_linalg_import)

        factor = phylogenetic_glm_module.cho_factor(
            matrix,
            lower=True,
            check_finite=False,
        )
        solved = phylogenetic_glm_module.cho_solve(
            factor,
            rhs,
            check_finite=False,
        )
    finally:
        phylogenetic_glm_module._CHO_FACTOR = previous_cho_factor
        phylogenetic_glm_module._CHO_SOLVE = previous_cho_solve

    np.testing.assert_allclose(matrix @ solved, rhs, atol=1e-12)


def test_normal_two_tailed_p_values_match_expected_values():
    p_values = phylogenetic_glm_module._normal_two_tailed_p_values(
        np.array([0.0, 1.959963984540054, -2.5758293035489004])
    )
    assert p_values == pytest.approx([1.0, 0.05, 0.01])


def test_normal_two_tailed_p_values_do_not_import_scipy_stats(monkeypatch):
    original_import = __import__

    def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "scipy.stats" or name.startswith("scipy.stats."):
            raise AssertionError("phylogenetic GLM z-test p-values should not import scipy.stats")
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr("builtins.__import__", fake_import)

    p_values = phylogenetic_glm_module._normal_two_tailed_p_values(
        np.array([0.0, 1.959963984540054])
    )
    assert p_values == pytest.approx([1.0, 0.05])


def test_normal_two_tailed_p_values_uses_vectorized_special_erfc(monkeypatch):
    captured = {}

    def fake_special_erfc(values):
        captured["values"] = values
        return np.full(values.shape, 0.25)

    monkeypatch.setattr(phylogenetic_glm_module, "special_erfc", fake_special_erfc)

    p_values = phylogenetic_glm_module._normal_two_tailed_p_values(
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
    A = rng.normal(size=(12, 12))
    info_matrix = A @ A.T + np.eye(12)
    scale = 1.7

    observed = phylogenetic_glm_module._standard_errors_from_info_matrix(
        info_matrix,
        scale=scale,
    )
    expected = np.sqrt(np.abs(np.diag(scale * np.linalg.inv(info_matrix))))

    np.testing.assert_allclose(observed, expected)


def test_standard_errors_from_info_matrix_spd_avoids_inverse(monkeypatch):
    rng = np.random.default_rng(20260628)
    A = rng.normal(size=(10, 10))
    info_matrix = A @ A.T + np.eye(10)

    def fail_inverse(*_args, **_kwargs):
        raise AssertionError("SPD information matrices should use Cholesky")

    monkeypatch.setattr(phylogenetic_glm_module.np.linalg, "inv", fail_inverse)

    observed = phylogenetic_glm_module._standard_errors_from_info_matrix(
        info_matrix
    )

    assert np.all(np.isfinite(observed))


def test_standard_errors_from_info_matrix_uses_diagonal_view(monkeypatch):
    rng = np.random.default_rng(20260701)
    A = rng.normal(size=(9, 9))
    info_matrix = A @ A.T + np.eye(9)
    scale = 1.4
    expected = np.sqrt(np.abs(scale * np.linalg.inv(info_matrix).diagonal()))

    def fail_diag(*_args, **_kwargs):
        raise AssertionError("standard errors should use ndarray diagonal access")

    monkeypatch.setattr(phylogenetic_glm_module.np, "diag", fail_diag)

    observed = phylogenetic_glm_module._standard_errors_from_info_matrix(
        info_matrix,
        scale=scale,
    )

    np.testing.assert_allclose(observed, expected)


def test_standard_errors_from_info_matrix_keeps_inverse_fallback():
    info_matrix = np.array([[1.0, 2.0], [2.0, 1.0]])

    observed = phylogenetic_glm_module._standard_errors_from_info_matrix(
        info_matrix
    )
    expected = np.sqrt(np.abs(np.diag(np.linalg.inv(info_matrix))))

    np.testing.assert_allclose(observed, expected)


@pytest.fixture
def binomial_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=GLM_TRAITS_FILE,
        response="binary_trait",
        predictors=["body_mass"],
        family="binomial",
        method=None,
        btol=10,
        log_alpha_bound=4,
        json=False,
    )


@pytest.fixture
def poisson_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=GLM_TRAITS_FILE,
        response="count_trait",
        predictors=["body_mass"],
        family="poisson",
        method=None,
        btol=10,
        log_alpha_bound=4,
        json=False,
    )


@pytest.fixture
def binomial_json_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=GLM_TRAITS_FILE,
        response="binary_trait",
        predictors=["body_mass"],
        family="binomial",
        method=None,
        btol=10,
        log_alpha_bound=4,
        json=True,
    )


@pytest.fixture
def poisson_json_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=GLM_TRAITS_FILE,
        response="count_trait",
        predictors=["body_mass"],
        family="poisson",
        method=None,
        btol=10,
        log_alpha_bound=4,
        json=True,
    )


@pytest.fixture
def multi_pred_binomial_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=GLM_TRAITS_FILE,
        response="binary_trait",
        predictors=["body_mass", "count_trait"],
        family="binomial",
        method=None,
        btol=10,
        log_alpha_bound=4,
        json=False,
    )


class TestProcessArgs:
    def test_defaults(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            response="y",
            predictors=["x1"],
            family="binomial",
        )
        svc = PhylogeneticGLM.__new__(PhylogeneticGLM)
        parsed = svc.process_args(args)
        assert parsed["tree_file_path"] == "t.tre"
        assert parsed["trait_data_path"] == "d.tsv"
        assert parsed["response"] == "y"
        assert parsed["predictors"] == ["x1"]
        assert parsed["family"] == "binomial"
        assert parsed["method"] == "logistic_MPLE"
        assert parsed["btol"] == 10
        assert parsed["log_alpha_bound"] == 4
        assert parsed["json_output"] is False

    def test_family_auto_method_binomial(self):
        args = Namespace(
            tree="t.tre", trait_data="d.tsv", response="y",
            predictors=["x1"], family="binomial", method=None,
        )
        svc = PhylogeneticGLM.__new__(PhylogeneticGLM)
        parsed = svc.process_args(args)
        assert parsed["method"] == "logistic_MPLE"

    def test_family_auto_method_poisson(self):
        args = Namespace(
            tree="t.tre", trait_data="d.tsv", response="y",
            predictors=["x1"], family="poisson", method=None,
        )
        svc = PhylogeneticGLM.__new__(PhylogeneticGLM)
        parsed = svc.process_args(args)
        assert parsed["method"] == "poisson_GEE"

    def test_custom_btol_and_log_alpha_bound(self):
        args = Namespace(
            tree="t.tre", trait_data="d.tsv", response="y",
            predictors=["x1"], family="binomial", method=None,
            btol=20, log_alpha_bound=8, json=True,
        )
        svc = PhylogeneticGLM.__new__(PhylogeneticGLM)
        parsed = svc.process_args(args)
        assert parsed["btol"] == 20
        assert parsed["log_alpha_bound"] == 8
        assert parsed["json_output"] is True


class TestValidation:
    def test_missing_column(self, binomial_args):
        binomial_args.response = "nonexistent_column"
        svc = PhylogeneticGLM(binomial_args)
        with pytest.raises(PhykitUserError) as exc_info:
            svc.run()
        assert any("nonexistent_column" in m for m in exc_info.value.messages)

    def test_response_in_predictors(self, binomial_args):
        binomial_args.predictors = ["body_mass", "binary_trait"]
        svc = PhylogeneticGLM(binomial_args)
        with pytest.raises(PhykitUserError) as exc_info:
            svc.run()
        assert any("must not also be a predictor" in m for m in exc_info.value.messages)

    def test_non_binary_for_binomial(self, poisson_args):
        poisson_args.family = "binomial"
        poisson_args.method = None
        poisson_args.response = "count_trait"
        svc = PhylogeneticGLM(poisson_args)
        with pytest.raises(PhykitUserError) as exc_info:
            svc.run()
        assert any("must contain only 0 and 1" in m for m in exc_info.value.messages)

    def test_insufficient_taxa(self, binomial_args):
        # Covered by _validate_tree requiring >= 3 tips
        pass

    def test_missing_predictor_column(self, binomial_args):
        binomial_args.predictors = ["nonexistent"]
        svc = PhylogeneticGLM(binomial_args)
        with pytest.raises(PhykitUserError) as exc_info:
            svc.run()
        assert any("nonexistent" in m for m in exc_info.value.messages)


class TestLogisticMPLE:
    def test_coefficients_reasonable(self, binomial_args):
        """Logistic MPLE should produce finite coefficients within btol bounds."""
        svc = PhylogeneticGLM(binomial_args)
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

        result = svc._fit_logistic_mple(tree, y, X, ordered_names)

        assert "coefficients" in result
        intercept = result["coefficients"]["(Intercept)"]["estimate"]
        slope = result["coefficients"]["body_mass"]["estimate"]
        assert np.isfinite(intercept)
        assert np.isfinite(slope)
        assert abs(intercept) <= svc.btol + 0.01
        assert abs(slope) <= svc.btol + 0.01

    def test_alpha_within_bounds(self, binomial_args):
        """Alpha should be within the bounds set by log_alpha_bound."""
        svc = PhylogeneticGLM(binomial_args)
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

        result = svc._fit_logistic_mple(tree, y, X, ordered_names)
        alpha = result["alpha"]
        assert alpha > 0
        assert alpha <= np.exp(svc.log_alpha_bound) + 0.01

    def test_standard_errors_positive(self, binomial_args):
        """Standard errors should be positive and finite."""
        svc = PhylogeneticGLM(binomial_args)
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

        result = svc._fit_logistic_mple(tree, y, X, ordered_names)
        for name, coef in result["coefficients"].items():
            assert coef["std_error"] > 0, f"SE for {name} should be positive"
            assert np.isfinite(coef["std_error"]), f"SE for {name} should be finite"

    def test_ultrametric_corrections(self, binomial_args):
        """_make_ultrametric should give non-negative D values."""
        svc = PhylogeneticGLM(binomial_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(GLM_TRAITS_FILE, tree_tips)
        ordered_names = sorted(traits.keys())

        D, Tmax, mean_height = svc._make_ultrametric(tree, ordered_names)
        assert Tmax > 0
        assert mean_height > 0
        assert mean_height <= Tmax
        assert np.all(D >= 0)
        # At least one tip should be at maximum distance (D=0)
        assert np.min(D) == pytest.approx(0.0)

    def test_root_tip_distances_fast_path_does_not_call_distance(
        self, binomial_args, monkeypatch
    ):
        svc = PhylogeneticGLM(binomial_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,C:3):0.5;"), "newick")

        def fail_depths(*args, **kwargs):
            raise AssertionError("depths fallback should not be called")

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("terminal fallback should not be called")

        def fail_distance(*args, **kwargs):
            raise AssertionError("distance fallback should not be called")

        def fail_array(*args, **kwargs):
            raise AssertionError("small standard-tree distances should use fromiter")

        monkeypatch.setattr(tree, "depths", fail_depths)
        monkeypatch.setattr(tree, "get_terminals", fail_get_terminals)
        monkeypatch.setattr(tree, "distance", fail_distance)
        monkeypatch.setattr(phylogenetic_glm_module.np, "array", fail_array)
        distances = svc._root_tip_distances(tree, ["A", "B", "C"])

        np.testing.assert_allclose(distances, np.array([2.0, 2.0, 3.0]))

    def test_large_root_tip_distances_avoid_reversed_iterator(
        self, binomial_args, monkeypatch
    ):
        from Bio.Phylo.BaseTree import Clade, Tree

        class NoReversedList(list):
            def __reversed__(self):
                raise AssertionError("large root-tip traversal should not call reversed")

        left = Clade(
            branch_length=1.0,
            clades=NoReversedList([
                Clade(name="A", branch_length=1.0),
                Clade(name="B", branch_length=1.0),
            ]),
        )
        root = Clade(
            clades=NoReversedList([left, Clade(name="C", branch_length=3.0)]),
        )
        tree = Tree(root=root)
        monkeypatch.setattr(
            phylogenetic_glm_module,
            "_ROOT_DISTANCE_BINARY_PUSH_MIN_TIPS",
            0,
        )

        distances = PhylogeneticGLM(binomial_args)._root_tip_distances(
            tree,
            ["A", "B", "C"],
        )

        np.testing.assert_allclose(distances, np.array([2.0, 2.0, 3.0]))

    def test_make_ultrametric_fast_path_does_not_call_distance(
        self, binomial_args, monkeypatch
    ):
        svc = PhylogeneticGLM(binomial_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,C:3):0.5;"), "newick")

        def fail_distance(*args, **kwargs):
            raise AssertionError("distance fallback should not be called")

        def fail_max(*args, **kwargs):
            raise AssertionError("ultrametric setup should use ndarray.max")

        def fail_mean(*args, **kwargs):
            raise AssertionError("ultrametric setup should use _mean_1d")

        monkeypatch.setattr(tree, "distance", fail_distance)
        monkeypatch.setattr(phylogenetic_glm_module.np, "max", fail_max)
        monkeypatch.setattr(phylogenetic_glm_module.np, "mean", fail_mean)
        D, Tmax, mean_height = svc._make_ultrametric(tree, ["A", "B", "C"])

        np.testing.assert_allclose(D, np.array([1.0, 1.0, 0.0]))
        assert Tmax == pytest.approx(3.0)
        assert mean_height == pytest.approx(7.0 / 3.0)

    def test_bernoulli_likelihood_finite(self, binomial_args):
        """Bernoulli log-likelihood should return a finite negative value."""
        svc = PhylogeneticGLM(binomial_args)
        tree = svc.read_tree_file()
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

        beta = np.array([0.0, 1.0])
        eta = X @ beta
        mu = 1.0 / (1.0 + np.exp(-eta))

        ll = PhylogeneticGLM._bernoulli_log_likelihood(y, mu)
        assert np.isfinite(ll)
        assert ll < 0  # Log-likelihood should be negative

    def test_bernoulli_likelihood_uses_array_sum(self, monkeypatch):
        y = np.array([0.0, 1.0, 1.0, 0.0])
        mu = np.array([0.2, 0.7, 0.8, 0.1])
        expected = float(np.sum(y * np.log(mu) + (1 - y) * np.log(1 - mu)))

        def fail_sum(*_args, **_kwargs):
            raise AssertionError("Bernoulli likelihood should use ndarray.sum")

        monkeypatch.setattr(phylogenetic_glm_module.np, "sum", fail_sum)

        assert PhylogeneticGLM._bernoulli_log_likelihood(y, mu) == pytest.approx(
            expected
        )

    def test_compute_dia_matches_scalar_reference_without_zero_loop(self, monkeypatch):
        svc = PhylogeneticGLM.__new__(PhylogeneticGLM)
        mu = np.array([0.12, 0.31, 0.47, 0.78, 0.93], dtype=float)
        D = np.array([0.0, 0.5, 1.25, 2.0, 3.5], dtype=float)
        alpha = 0.17
        meanp = float(np.mean(mu))
        meanq = 1.0 - meanp
        expected = []
        for mui, Di in zip(mu, D):
            if mui < meanp:
                m = mui * np.sqrt(meanq / max(meanp, 1e-300))
            else:
                m = (1.0 - mui) * np.sqrt(meanp / max(meanq, 1e-300))
            expected.append(np.sqrt(max(m * m, 1e-300)) * np.exp(alpha * Di))

        def fail_zeros(*_args, **_kwargs):
            raise AssertionError("dia computation should use vectorized selection")

        monkeypatch.setattr(phylogenetic_glm_module.np, "zeros", fail_zeros)

        observed = svc._compute_dia(["A", "B", "C", "D", "E"], alpha, mu, D)

        np.testing.assert_allclose(observed, np.array(expected))

    def test_logistic_starting_values_match_diagonal_weight_reference(self, monkeypatch):
        svc = PhylogeneticGLM.__new__(PhylogeneticGLM)
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

        expected = reference()

        def fail_max(*_args, **_kwargs):
            raise AssertionError("starting-value convergence should use ndarray.max")

        monkeypatch.setattr(phylogenetic_glm_module.np, "max", fail_max)

        observed = svc._logistic_starting_values(y, X, btol=10)
        np.testing.assert_allclose(observed, expected)

class TestPoissonGEE:
    def test_convergence(self, poisson_args):
        """Poisson GEE should converge for our test data."""
        svc = PhylogeneticGLM(poisson_args)
        tree = svc.read_tree_file()
        svc.validate_tree(tree, min_tips=3, require_branch_lengths=True)
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(GLM_TRAITS_FILE, tree_tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        resp_idx = trait_names.index("count_trait")
        pred_idx = trait_names.index("body_mass")
        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X = np.column_stack([
            np.ones(n),
            np.array([traits[name][pred_idx] for name in ordered_names]),
        ])

        result = svc._fit_poisson_gee(tree, y, X, ordered_names)

        # Should have sensible coefficients
        intercept = result["coefficients"]["(Intercept)"]["estimate"]
        slope = result["coefficients"]["body_mass"]["estimate"]
        assert np.isfinite(intercept)
        assert np.isfinite(slope)

    def test_positive_slope(self, poisson_args):
        """Count trait increases with body_mass, so slope should be positive."""
        svc = PhylogeneticGLM(poisson_args)
        tree = svc.read_tree_file()
        svc.validate_tree(tree, min_tips=3, require_branch_lengths=True)
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(GLM_TRAITS_FILE, tree_tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        resp_idx = trait_names.index("count_trait")
        pred_idx = trait_names.index("body_mass")
        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X = np.column_stack([
            np.ones(n),
            np.array([traits[name][pred_idx] for name in ordered_names]),
        ])

        result = svc._fit_poisson_gee(tree, y, X, ordered_names)
        slope = result["coefficients"]["body_mass"]["estimate"]
        assert slope > 0

    def test_overdispersion_positive(self, poisson_args):
        """Overdispersion phi should be positive."""
        svc = PhylogeneticGLM(poisson_args)
        tree = svc.read_tree_file()
        svc.validate_tree(tree, min_tips=3, require_branch_lengths=True)
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(GLM_TRAITS_FILE, tree_tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        resp_idx = trait_names.index("count_trait")
        pred_idx = trait_names.index("body_mass")
        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X = np.column_stack([
            np.ones(n),
            np.array([traits[name][pred_idx] for name in ordered_names]),
        ])

        result = svc._fit_poisson_gee(tree, y, X, ordered_names)
        assert result["overdispersion"] > 0

    def test_standard_errors_positive(self, poisson_args):
        """Standard errors should be positive."""
        svc = PhylogeneticGLM(poisson_args)
        tree = svc.read_tree_file()
        svc.validate_tree(tree, min_tips=3, require_branch_lengths=True)
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(GLM_TRAITS_FILE, tree_tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        resp_idx = trait_names.index("count_trait")
        pred_idx = trait_names.index("body_mass")
        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X = np.column_stack([
            np.ones(n),
            np.array([traits[name][pred_idx] for name in ordered_names]),
        ])

        result = svc._fit_poisson_gee(tree, y, X, ordered_names)
        for name, coef in result["coefficients"].items():
            assert coef["std_error"] > 0, f"SE for {name} should be positive"

    def test_starting_values(self, poisson_args):
        """Starting values should give reasonable initial estimates."""
        svc = PhylogeneticGLM(poisson_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(GLM_TRAITS_FILE, tree_tips)
        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        resp_idx = trait_names.index("count_trait")
        pred_idx = trait_names.index("body_mass")
        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X = np.column_stack([
            np.ones(n),
            np.array([traits[name][pred_idx] for name in ordered_names]),
        ])

        beta0 = svc._poisson_starting_values(y, X)
        assert np.all(np.isfinite(beta0))
        # Intercept should be near log(mean(y))
        assert abs(beta0[0] - np.log(np.mean(y))) < 2.0

    def test_poisson_starting_values_match_diagonal_weight_reference(self, monkeypatch):
        svc = PhylogeneticGLM.__new__(PhylogeneticGLM)
        rng = np.random.default_rng(20260623)
        n = 40
        X = np.column_stack([np.ones(n), rng.normal(size=(n, 3))])
        y = rng.poisson(2.5, size=n).astype(float)

        def reference():
            beta = np.zeros(X.shape[1])
            beta[0] = np.log(max(np.mean(y), 0.1))
            for _ in range(50):
                eta = np.clip(X @ beta, -20, 20)
                mu = np.exp(eta)
                mu = np.clip(mu, 1e-6, 1e10)
                weights = np.diag(mu)
                z = eta + (y - mu) / mu
                beta_new = np.linalg.solve(X.T @ weights @ X, X.T @ weights @ z)
                if np.max(np.abs(beta_new - beta)) < 1e-8:
                    return beta_new
                beta = beta_new
            return beta

        expected = reference()

        def fail_max(*_args, **_kwargs):
            raise AssertionError("starting-value convergence should use ndarray.max")

        monkeypatch.setattr(phylogenetic_glm_module.np, "max", fail_max)

        observed = svc._poisson_starting_values(y, X)
        np.testing.assert_allclose(observed, expected)

    def test_poisson_gee_information_and_score_match_diagonal_reference(self):
        rng = np.random.default_rng(20260623)
        n = 12
        X = np.column_stack([np.ones(n), rng.normal(size=(n, 2))])
        A = rng.normal(size=(n, n))
        R = A @ A.T + np.eye(n)
        R_inv = np.linalg.inv(R)
        mu = np.exp(rng.normal(size=n) * 0.2)
        y = rng.poisson(mu).astype(float)

        I_mat, score = PhylogeneticGLM._poisson_gee_information_and_score(
            X, R_inv, mu, y
        )

        sqrt_mu = np.sqrt(mu)
        WX = np.diag(sqrt_mu) @ X
        expected_I = WX.T @ R_inv @ WX
        rhs = (y - mu) / sqrt_mu
        expected_score = WX.T @ (R_inv @ rhs)

        np.testing.assert_allclose(I_mat, expected_I)
        np.testing.assert_allclose(score, expected_score)

    def test_poisson_gee_information_and_score_cholesky_matches_inverse(self):
        rng = np.random.default_rng(20260628)
        n = 18
        X = np.column_stack([np.ones(n), rng.normal(size=(n, 3))])
        A = rng.normal(size=(n, n))
        R = A @ A.T + np.eye(n)
        R_inv = np.linalg.inv(R)
        R_factor = phylogenetic_glm_module.cho_factor(
            R, lower=True, check_finite=False
        )
        mu = np.exp(rng.normal(size=n) * 0.2)
        y = rng.poisson(mu).astype(float)

        expected_I, expected_score = PhylogeneticGLM._poisson_gee_information_and_score(
            X, R_inv, mu, y
        )
        observed_I, observed_score = (
            PhylogeneticGLM._poisson_gee_information_and_score_cholesky(
                X, R_factor, mu, y
            )
        )

        np.testing.assert_allclose(observed_I, expected_I)
        np.testing.assert_allclose(observed_score, expected_score)

    def test_poisson_gee_information_cholesky_matches_score_helper(self):
        rng = np.random.default_rng(20260628)
        n = 18
        X = np.column_stack([np.ones(n), rng.normal(size=(n, 3))])
        A = rng.normal(size=(n, n))
        R = A @ A.T + np.eye(n)
        R_factor = phylogenetic_glm_module.cho_factor(
            R, lower=True, check_finite=False
        )
        mu = np.exp(rng.normal(size=n) * 0.2)

        expected_I, _ = PhylogeneticGLM._poisson_gee_information_and_score_cholesky(
            X, R_factor, mu, mu
        )
        observed_I = PhylogeneticGLM._poisson_gee_information_cholesky(
            X, R_factor, mu
        )

        np.testing.assert_allclose(observed_I, expected_I)

    def test_poisson_gee_information_cholesky_solves_design_columns_only(
        self, monkeypatch
    ):
        rng = np.random.default_rng(20260628)
        n = 18
        X = np.column_stack([np.ones(n), rng.normal(size=(n, 3))])
        A = rng.normal(size=(n, n))
        R = A @ A.T + np.eye(n)
        R_factor = phylogenetic_glm_module.cho_factor(
            R, lower=True, check_finite=False
        )
        mu = np.exp(rng.normal(size=n) * 0.2)
        original_cho_solve = phylogenetic_glm_module.cho_solve
        rhs_shapes = []

        def recording_cho_solve(factor, rhs, **kwargs):
            rhs_shapes.append(rhs.shape)
            return original_cho_solve(factor, rhs, **kwargs)

        monkeypatch.setattr(phylogenetic_glm_module, "cho_solve", recording_cho_solve)

        PhylogeneticGLM._poisson_gee_information_cholesky(X, R_factor, mu)

        assert rhs_shapes == [(n, X.shape[1])]

    def test_poisson_gee_information_and_score_cholesky_uses_single_solve(
        self, monkeypatch
    ):
        rng = np.random.default_rng(20260628)
        n = 18
        X = np.column_stack([np.ones(n), rng.normal(size=(n, 3))])
        A = rng.normal(size=(n, n))
        R = A @ A.T + np.eye(n)
        R_factor = phylogenetic_glm_module.cho_factor(
            R, lower=True, check_finite=False
        )
        mu = np.exp(rng.normal(size=n) * 0.2)
        y = rng.poisson(mu).astype(float)
        original_cho_solve = phylogenetic_glm_module.cho_solve
        solve_calls = 0

        def counting_cho_solve(*args, **kwargs):
            nonlocal solve_calls
            solve_calls += 1
            return original_cho_solve(*args, **kwargs)

        monkeypatch.setattr(phylogenetic_glm_module, "cho_solve", counting_cho_solve)

        PhylogeneticGLM._poisson_gee_information_and_score_cholesky(
            X, R_factor, mu, y
        )

        assert solve_calls == 1

    def test_fit_poisson_gee_uses_cholesky_information_path(
        self, poisson_args, monkeypatch
    ):
        svc = PhylogeneticGLM(poisson_args)
        ordered_names = ["A", "B", "C", "D"]
        y = np.array([1.0, 2.0, 3.0, 4.0])
        X = np.column_stack([np.ones(4), np.array([0.1, 0.2, 0.3, 0.4])])
        svc._precomputed_vcv = np.array(
            [
                [1.0, 0.2, 0.1, 0.1],
                [0.2, 1.0, 0.1, 0.1],
                [0.1, 0.1, 1.0, 0.2],
                [0.1, 0.1, 0.2, 1.0],
            ]
        )

        def fail_inverse_information(*_args, **_kwargs):
            raise AssertionError("positive-definite R should use Cholesky path")

        def fail_diag(*_args, **_kwargs):
            raise AssertionError("Poisson GEE should use ndarray diagonal access")

        monkeypatch.setattr(
            svc,
            "_poisson_gee_information_and_score",
            fail_inverse_information,
        )
        monkeypatch.setattr(phylogenetic_glm_module.np, "diag", fail_diag)

        result = svc._fit_poisson_gee(tree=None, y=y, X=X, ordered_names=ordered_names)

        assert result["family"] == "poisson"
        assert all(np.isfinite(row["estimate"]) for row in result["coefficients"].values())

    def test_fit_poisson_gee_convergence_uses_delta_array_sum(
        self, poisson_args, monkeypatch
    ):
        svc = PhylogeneticGLM(poisson_args)
        ordered_names = ["A", "B", "C", "D"]
        y = np.array([1.0, 2.0, 3.0, 4.0])
        X = np.column_stack([np.ones(4), np.array([0.1, 0.2, 0.3, 0.4])])
        svc._precomputed_vcv = np.eye(4)

        def tiny_step(*_args, **_kwargs):
            return np.eye(2), np.array([1e-12, -2e-12])

        def fail_delta_sum(values, *args, **kwargs):
            arr = np.asarray(values)
            if arr.shape == (2,):
                raise AssertionError("Poisson GEE convergence should use delta.sum")
            return original_sum(values, *args, **kwargs)

        original_sum = phylogenetic_glm_module.np.sum
        monkeypatch.setattr(
            svc,
            "_poisson_gee_information_and_score_cholesky",
            tiny_step,
        )
        monkeypatch.setattr(
            svc,
            "_poisson_gee_information_cholesky",
            lambda *_args, **_kwargs: np.eye(2),
        )
        monkeypatch.setattr(phylogenetic_glm_module.np, "sum", fail_delta_sum)

        result = svc._fit_poisson_gee(tree=None, y=y, X=X, ordered_names=ordered_names)

        assert result["family"] == "poisson"

    def test_fit_poisson_gee_uses_information_only_final_covariance(
        self, poisson_args, monkeypatch
    ):
        svc = PhylogeneticGLM(poisson_args)
        ordered_names = ["A", "B", "C", "D"]
        y = np.array([1.0, 2.0, 3.0, 4.0])
        X = np.column_stack([np.ones(4), np.array([0.1, 0.2, 0.3, 0.4])])
        svc._precomputed_vcv = np.array(
            [
                [1.0, 0.2, 0.1, 0.1],
                [0.2, 1.0, 0.1, 0.1],
                [0.1, 0.1, 1.0, 0.2],
                [0.1, 0.1, 0.2, 1.0],
            ]
        )
        original_information = svc._poisson_gee_information_cholesky
        final_information_calls = 0

        def counting_information(*args, **kwargs):
            nonlocal final_information_calls
            final_information_calls += 1
            return original_information(*args, **kwargs)

        monkeypatch.setattr(
            svc,
            "_poisson_gee_information_cholesky",
            counting_information,
        )

        result = svc._fit_poisson_gee(tree=None, y=y, X=X, ordered_names=ordered_names)

        assert result["family"] == "poisson"
        assert final_information_calls == 1

class TestRun:
    def test_run_uses_unmodified_tree_read(self, mocker):
        args = Namespace(
            tree="/some/path/to/file.tre",
            trait_data="/some/path/to/traits.tsv",
            response="y",
            predictors=["x"],
            family="binomial",
            method="logistic_MPLE",
            json=False,
            gene_trees=None,
            btol=10,
            log_alpha_bound=4,
        )
        svc = PhylogeneticGLM(args)
        tree = object()
        mocked_read = mocker.patch.object(
            PhylogeneticGLM,
            "read_tree_file_unmodified",
            return_value=tree,
        )
        mocker.patch.object(
            PhylogeneticGLM,
            "read_tree_file",
            side_effect=AssertionError("run should not copy cached trees"),
        )
        mocked_validate = mocker.patch.object(
            PhylogeneticGLM, "validate_tree"
        )
        mocker.patch.object(
            PhylogeneticGLM,
            "get_tip_names_from_tree",
            return_value=["a", "b", "c"],
        )
        mocker.patch(
            "phykit.services.tree.phylogenetic_glm.parse_multi_trait_file",
            return_value=(
                ["y", "x"],
                {
                    "a": [0.0, 1.0],
                    "b": [1.0, 2.0],
                    "c": [0.0, 3.0],
                },
            ),
        )
        fit_result = {"log_likelihood": -1.0}
        null_result = {"log_likelihood": -1.5}
        mocked_fit = mocker.patch.object(
            PhylogeneticGLM,
            "_fit_logistic_mple",
            side_effect=[fit_result, null_result],
        )
        mocker.patch.object(PhylogeneticGLM, "_print_text_output")

        svc.run()

        mocked_read.assert_called_once_with()
        mocked_validate.assert_called_once_with(
            tree,
            min_tips=3,
            require_branch_lengths=True,
            context="phylogenetic GLM",
        )
        assert mocked_fit.call_args_list[0].args[0] is tree
        assert mocked_fit.call_args_list[1].args[0] is tree

    def test_run_resolves_trait_columns_with_first_index_map(self, mocker):
        args = Namespace(
            tree="/some/path/to/file.tre",
            trait_data="/some/path/to/traits.tsv",
            response="y",
            predictors=["x"],
            family="binomial",
            method="logistic_MPLE",
            json=False,
            gene_trees=None,
            btol=10,
            log_alpha_bound=4,
        )
        svc = PhylogeneticGLM(args)
        tree = object()
        mocker.patch.object(svc, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(svc, "validate_tree")
        mocker.patch.object(svc, "get_tip_names_from_tree", return_value=["a", "b", "c"])
        mocker.patch(
            "phykit.services.tree.phylogenetic_glm.parse_multi_trait_file",
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
            "phykit.services.tree.phylogenetic_glm.response_predictor_arrays",
            side_effect=arrays,
        )
        mocker.patch.object(
            PhylogeneticGLM,
            "_fit_logistic_mple",
            side_effect=[{"log_likelihood": -1.0}, {"log_likelihood": -1.5}],
        )
        mocker.patch.object(PhylogeneticGLM, "_print_text_output")

        svc.run()

    @patch("builtins.print")
    def test_text_output_binomial(self, mocked_print, binomial_args):
        svc = PhylogeneticGLM(binomial_args)
        svc.run()
        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Phylogenetic GLM (Logistic MPLE)" in all_output
        assert "(Intercept)" in all_output
        assert "body_mass" in all_output
        assert "Estimated alpha" in all_output
        assert "binary_trait ~ body_mass" in all_output

    @patch("builtins.print")
    def test_text_output_poisson(self, mocked_print, poisson_args):
        svc = PhylogeneticGLM(poisson_args)
        svc.run()
        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Phylogenetic GLM (Poisson GEE)" in all_output
        assert "(Intercept)" in all_output
        assert "body_mass" in all_output
        assert "Overdispersion (phi)" in all_output
        assert "count_trait ~ body_mass" in all_output

    @patch("builtins.print")
    def test_json_output_binomial(self, mocked_print, binomial_json_args):
        svc = PhylogeneticGLM(binomial_json_args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["family"] == "binomial"
        assert payload["method"] == "logistic_MPLE"
        assert "alpha" in payload
        assert "(Intercept)" in payload["coefficients"]
        assert "body_mass" in payload["coefficients"]
        assert payload["formula"] == "binary_trait ~ body_mass"

    @patch("builtins.print")
    def test_json_output_poisson(self, mocked_print, poisson_json_args):
        svc = PhylogeneticGLM(poisson_json_args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["family"] == "poisson"
        assert payload["method"] == "poisson_GEE"
        assert "overdispersion" in payload
        assert "(Intercept)" in payload["coefficients"]
        assert payload["formula"] == "count_trait ~ body_mass"

    def test_format_result_maps_parallel_arrays_by_name(self):
        svc = PhylogeneticGLM.__new__(PhylogeneticGLM)

        result = svc._format_result(
            family="poisson",
            method="poisson_GEE",
            coef_names=["(Intercept)", "body_mass", "count_trait"],
            beta_hat=np.array([0.25, 1.5, -0.75]),
            se=np.array([0.1, 0.2, 0.3]),
            z_stats=np.array([2.5, 7.5, -2.5]),
            p_values=np.array([0.03, 0.001, 0.04]),
            ll=-3.5,
            aic=15.0,
            formula="count_trait ~ body_mass + count_trait",
            n=3,
            k=2,
            ordered_names=["taxon_a", "taxon_b", "taxon_c"],
            fitted={
                "taxon_a": 1.1,
                "taxon_b": 2.2,
                "taxon_c": 3.3,
            },
            alpha=0.5,
            overdispersion=1.25,
        )

        assert result["coefficients"]["body_mass"] == {
            "estimate": 1.5,
            "std_error": 0.2,
            "z_value": 7.5,
            "p_value": 0.001,
        }
        assert result["fitted_values"] == {
            "taxon_a": 1.1,
            "taxon_b": 2.2,
            "taxon_c": 3.3,
        }
        assert result["overdispersion"] == 1.25

    @patch("builtins.print")
    def test_signif_codes_displayed(self, mocked_print, poisson_args):
        svc = PhylogeneticGLM(poisson_args)
        svc.run()
        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Signif. codes:" in all_output

    @patch("builtins.print")
    def test_print_text_output_batches_coefficient_rows(self, mocked_print):
        svc = PhylogeneticGLM.__new__(PhylogeneticGLM)
        result = {
            "family": "binomial",
            "method": "logistic_MPLE",
            "formula": "binary_trait ~ body_mass + count_trait",
            "alpha": 0.125,
            "overdispersion": None,
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
            "aic": 31.0,
            "pseudo_r_squared_mcfadden": 0.42,
            "n_observations": 25,
        }

        svc._print_text_output(result)

        mocked_print.assert_called_once()
        output = mocked_print.call_args.args[0]
        assert "Phylogenetic GLM (Logistic MPLE)" in output
        assert "binary_trait ~ body_mass + count_trait" in output
        assert "(Intercept)" in output
        assert "body_mass" in output
        assert "    ***" in output
        assert "    *" in output
        assert "Pseudo-R² (McFadden): 0.4200" in output
        assert "Number of observations: 25" in output

    @patch("builtins.print")
    def test_multiple_predictors(self, mocked_print, multi_pred_binomial_args):
        svc = PhylogeneticGLM(multi_pred_binomial_args)
        svc.run()
        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "body_mass" in all_output
        assert "count_trait" in all_output
        assert "binary_trait ~ body_mass + count_trait" in all_output


GENE_TREES_FILE = str(SAMPLE_FILES / "gene_trees_simple.nwk")


class TestDiscordanceVCV:
    def test_run_reuses_single_discordance_vcv_build(self, mocker):
        args = Namespace(
            tree="tree.tre",
            trait_data="traits.tsv",
            response="count_trait",
            predictors=["body_mass"],
            family="poisson",
            method=None,
            btol=10,
            log_alpha_bound=4,
            json=False,
            gene_trees="gene_trees.nwk",
        )
        svc = PhylogeneticGLM(args)
        tree = object()
        vcv = np.eye(3)
        fit_ordered_names = []
        precomputed_vcvs = []

        mocker.patch.object(svc, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(svc, "validate_tree")
        mocker.patch.object(
            svc,
            "get_tip_names_from_tree",
            return_value=["a", "b", "c", "d"],
        )
        mocker.patch(
            "phykit.services.tree.phylogenetic_glm.parse_multi_trait_file",
            return_value=(
                ["count_trait", "body_mass"],
                {
                    "a": [1.0, 10.0],
                    "b": [2.0, 20.0],
                    "c": [3.0, 30.0],
                    "d": [4.0, 40.0],
                },
            ),
        )

        import phykit.services.tree.vcv_utils as vcv_utils

        parse_gene_trees = mocker.patch.object(
            vcv_utils,
            "parse_gene_trees",
            return_value=["g1", "g2"],
        )

        def build_discordance_vcv(tree_arg, gene_trees, ordered_names):
            assert tree_arg is tree
            assert gene_trees == ["g1", "g2"]
            assert ordered_names == ["a", "b", "c", "d"]
            return vcv, {
                "shared_taxa": ["a", "c", "d"],
                "n_gene_trees": 2,
                "n_shared_taxa": 3,
            }

        build_vcv = mocker.patch.object(
            vcv_utils,
            "build_discordance_vcv",
            side_effect=build_discordance_vcv,
        )

        def fit_poisson(_tree, _y, _X, ordered_names):
            fit_ordered_names.append(list(ordered_names))
            precomputed_vcvs.append(svc._precomputed_vcv)
            return {"log_likelihood": -1.0}

        mocker.patch.object(svc, "_fit_poisson_gee", side_effect=fit_poisson)
        mocker.patch.object(svc, "_print_text_output")

        svc.run()

        parse_gene_trees.assert_called_once_with("gene_trees.nwk")
        build_vcv.assert_called_once()
        assert fit_ordered_names == [["a", "c", "d"], ["a", "c", "d"]]
        assert precomputed_vcvs == [vcv, vcv]

    def test_poisson_gee_with_gene_trees_json(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=GLM_TRAITS_FILE,
            response="count_trait",
            predictors=["body_mass"],
            family="poisson",
            method=None,
            btol=10,
            log_alpha_bound=4,
            json=True,
            gene_trees=GENE_TREES_FILE,
        )
        svc = PhylogeneticGLM(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "vcv_metadata" in data
        assert data["vcv_metadata"]["n_gene_trees"] == 10
        assert "coefficients" in data

    def test_binomial_with_gene_trees_json(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=GLM_TRAITS_FILE,
            response="binary_trait",
            predictors=["body_mass"],
            family="binomial",
            method=None,
            btol=10,
            log_alpha_bound=4,
            json=True,
            gene_trees=GENE_TREES_FILE,
        )
        svc = PhylogeneticGLM(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "vcv_metadata" in data
        assert "coefficients" in data

    def test_discordance_changes_poisson_coefficients(self, capsys):
        """With discordant gene trees, Poisson GEE coefficients should differ."""
        base_kwargs = dict(
            tree=TREE_SIMPLE,
            trait_data=GLM_TRAITS_FILE,
            response="count_trait",
            predictors=["body_mass"],
            family="poisson",
            method=None,
            btol=10,
            log_alpha_bound=4,
            json=True,
        )
        args_no_gt = Namespace(**base_kwargs)
        svc = PhylogeneticGLM(args_no_gt)
        svc.run()
        out, _ = capsys.readouterr()
        coefs_no_gt = json.loads(out)["coefficients"]

        args_gt = Namespace(**base_kwargs, gene_trees=GENE_TREES_FILE)
        svc = PhylogeneticGLM(args_gt)
        svc.run()
        out, _ = capsys.readouterr()
        coefs_gt = json.loads(out)["coefficients"]

        # Coefficients should differ with discordant gene trees
        assert coefs_gt["body_mass"]["estimate"] != pytest.approx(
            coefs_no_gt["body_mass"]["estimate"], abs=1e-6
        )

    @patch("builtins.print")
    def test_text_output_with_gene_trees(self, mocked_print):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=GLM_TRAITS_FILE,
            response="count_trait",
            predictors=["body_mass"],
            family="poisson",
            method=None,
            btol=10,
            log_alpha_bound=4,
            json=False,
            gene_trees=GENE_TREES_FILE,
        )
        svc = PhylogeneticGLM(args)
        svc.run()
        assert mocked_print.called


class TestEffectSize:
    def test_pseudo_r2_poisson_json(self, capsys):
        """JSON output for Poisson should contain pseudo_r_squared_mcfadden and ll_null."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=GLM_TRAITS_FILE,
            response="count_trait",
            predictors=["body_mass"],
            family="poisson",
            method=None,
            btol=10,
            log_alpha_bound=4,
            json=True,
        )
        svc = PhylogeneticGLM(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "pseudo_r_squared_mcfadden" in data
        assert "ll_null" in data
        assert data["pseudo_r_squared_mcfadden"] >= 0

    def test_pseudo_r2_binomial_json(self, capsys):
        """JSON output for binomial should contain pseudo_r_squared_mcfadden and ll_null."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=GLM_TRAITS_FILE,
            response="binary_trait",
            predictors=["body_mass"],
            family="binomial",
            method=None,
            btol=10,
            log_alpha_bound=4,
            json=True,
        )
        svc = PhylogeneticGLM(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "pseudo_r_squared_mcfadden" in data
        assert "ll_null" in data
        assert data["pseudo_r_squared_mcfadden"] >= 0

    def test_pseudo_r2_between_0_and_1(self, capsys):
        """Pseudo-R² for Poisson should be between 0 and 1."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=GLM_TRAITS_FILE,
            response="count_trait",
            predictors=["body_mass"],
            family="poisson",
            method=None,
            btol=10,
            log_alpha_bound=4,
            json=True,
        )
        svc = PhylogeneticGLM(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        r2 = data["pseudo_r_squared_mcfadden"]
        assert 0 <= r2 <= 1, f"Pseudo-R² should be in [0, 1], got {r2}"

    @patch("builtins.print")
    def test_pseudo_r2_in_text_output(self, mocked_print):
        """Text output should include Pseudo-R² line."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=GLM_TRAITS_FILE,
            response="count_trait",
            predictors=["body_mass"],
            family="poisson",
            method=None,
            btol=10,
            log_alpha_bound=4,
            json=False,
        )
        svc = PhylogeneticGLM(args)
        svc.run()
        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Pseudo-R" in all_output

    def test_pseudo_r2_poisson_matches_reference(self, capsys):
        """Poisson McFadden's pseudo-R² must match reference values.

        Reference (PhyKIT Poisson GEE, cross-checked with formula
        1 - LL_full/LL_null):
          ll_full  = -13.4023801615
          ll_null  = -18.0173819720
          pseudo_r2 = 1 - (-13.4024 / -18.0174) = 0.2561416424
        """
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=GLM_TRAITS_FILE,
            response="count_trait",
            predictors=["body_mass"],
            family="poisson",
            method=None,
            btol=10,
            log_alpha_bound=4,
            json=True,
        )
        svc = PhylogeneticGLM(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert data["pseudo_r_squared_mcfadden"] == pytest.approx(
            0.2561416424, abs=1e-4
        )
        assert data["ll_null"] == pytest.approx(-18.0173819720, abs=1e-2)

    def test_pseudo_r2_binomial_matches_reference(self, capsys):
        """Binomial McFadden's pseudo-R² must match reference values.

        Reference (PhyKIT logistic MPLE):
          ll_full  = -1.8326609276
          ll_null  = -4.7365924439
          pseudo_r2 = 1 - (-1.8327 / -4.7366) = 0.6130845224
        """
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=GLM_TRAITS_FILE,
            response="binary_trait",
            predictors=["body_mass"],
            family="binomial",
            method=None,
            btol=10,
            log_alpha_bound=4,
            json=True,
        )
        svc = PhylogeneticGLM(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert data["pseudo_r_squared_mcfadden"] == pytest.approx(
            0.6130845224, abs=1e-4
        )
        assert data["ll_null"] == pytest.approx(-4.7365924439, abs=1e-2)
