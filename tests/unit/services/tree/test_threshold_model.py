import builtins
import importlib
import json
import os
import subprocess
import sys
from argparse import Namespace
from io import StringIO
from pathlib import Path

import numpy as np
import pytest
from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin

import phykit.services.tree.threshold_model as threshold_model_module
from phykit.services.tree.threshold_model import ThresholdModel
from phykit.errors import PhykitUserError


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_threshold_traits.tsv")


def test_module_import_does_not_import_scipy_special_or_linalg(monkeypatch):
    module_name = "phykit.services.tree.threshold_model"
    previous = sys.modules.pop(module_name, None)
    original_import = builtins.__import__

    def guarded_import(name, globals=None, locals=None, fromlist=(), level=0):
        if (
            name == "scipy.special"
            or name.startswith("scipy.special.")
            or name == "scipy.linalg"
            or name.startswith("scipy.linalg.")
        ):
            raise AssertionError(
                "threshold_model module import should not import SciPy special/linalg"
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


def test_module_import_does_not_import_numpy_or_scipy_special():
    code = """
import sys
import phykit.services.tree.threshold_model as module
assert hasattr(module.np, "__getattr__")
assert callable(module.print_json)
assert "typing" not in sys.modules
assert "pickle" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
assert "phykit.services.tree.vcv_utils" not in sys.modules
assert "numpy" not in sys.modules
assert "scipy.special" not in sys.modules
assert "scipy.linalg" not in sys.modules
assert "scipy.stats" not in sys.modules
assert module._CHO_FACTOR is None
assert module._CHO_SOLVE is None
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_vcv_inverse_and_logdet_cholesky_matches_inverse_and_slogdet():
    rng = np.random.default_rng(20260628)
    A = rng.normal(size=(8, 8))
    C = A @ A.T + np.eye(8) * 0.25

    observed_inv, observed_logdet = ThresholdModel._vcv_inverse_and_logdet(C)
    expected_inv = np.linalg.inv(C)
    _, expected_logdet = np.linalg.slogdet(C)

    np.testing.assert_allclose(observed_inv, expected_inv)
    assert observed_logdet == pytest.approx(expected_logdet)


def test_vcv_inverse_and_logdet_spd_avoids_inverse_and_slogdet(monkeypatch):
    rng = np.random.default_rng(20260628)
    A = rng.normal(size=(6, 6))
    C = A @ A.T + np.eye(6)

    def fail_inverse(_matrix):
        raise AssertionError("SPD VCV setup should use Cholesky")

    def fail_slogdet(_matrix):
        raise AssertionError("SPD VCV logdet should use the Cholesky factor")

    monkeypatch.setattr(threshold_model_module.np.linalg, "inv", fail_inverse)
    monkeypatch.setattr(threshold_model_module.np.linalg, "slogdet", fail_slogdet)

    observed_inv, observed_logdet = ThresholdModel._vcv_inverse_and_logdet(C)

    assert observed_inv.shape == C.shape
    assert np.all(np.isfinite(observed_inv))
    assert np.isfinite(observed_logdet)


def test_vcv_inverse_and_logdet_keeps_inverse_fallback():
    C = np.array([[-2.0, 0.0], [0.0, -1.0]])

    observed_inv, observed_logdet = ThresholdModel._vcv_inverse_and_logdet(C)
    expected_inv = np.linalg.inv(C)
    _, expected_logdet = np.linalg.slogdet(C)

    np.testing.assert_allclose(observed_inv, expected_inv)
    assert observed_logdet == pytest.approx(expected_logdet)


def _make_tree():
    return Phylo.read(
        StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick"
    )


class TestProcessArgs:
    def test_defaults(self):
        args = Namespace(
            tree="t.tre",
            trait_data="traits.tsv",
            traits="habitat,body_mass",
            types="discrete,continuous",
        )
        svc = ThresholdModel(args)
        assert svc.tree_file_path == "t.tre"
        assert svc.traits == ["habitat", "body_mass"]
        assert svc.types == ["discrete", "continuous"]
        assert svc.ngen == 100000
        assert svc.sample == 100
        assert svc.burnin == 0.2
        assert svc.seed is None
        assert svc.plot_output is None
        assert svc.json_output is False

    def test_all_options(self):
        args = Namespace(
            tree="t.tre",
            trait_data="traits.tsv",
            traits="a,b",
            types="continuous,continuous",
            ngen=5000,
            sample=10,
            burnin=0.1,
            seed=42,
            plot_output="trace.png",
            json=True,
        )
        svc = ThresholdModel(args)
        assert svc.ngen == 5000
        assert svc.sample == 10
        assert svc.burnin == 0.1
        assert svc.seed == 42
        assert svc.plot_output == "trace.png"
        assert svc.json_output is True

    def test_json_flag(self):
        args = Namespace(
            tree="t.tre",
            trait_data="traits.tsv",
            traits="a,b",
            types="discrete,discrete",
            json=True,
        )
        svc = ThresholdModel(args)
        assert svc.json_output is True

    def test_missing_traits_raises(self):
        args = Namespace(
            tree="t.tre",
            trait_data="traits.tsv",
            traits=None,
            types="discrete,continuous",
        )
        with pytest.raises(SystemExit):
            ThresholdModel(args)

    def test_invalid_type_raises(self):
        args = Namespace(
            tree="t.tre",
            trait_data="traits.tsv",
            traits="a,b",
            types="discrete,ordinal",
        )
        with pytest.raises(SystemExit):
            ThresholdModel(args)


class TestParseMultiTraitFile:
    def test_valid_parse(self, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "taxon\tt1\tt2\n"
            "A\t0\t1.5\n"
            "B\t1\t2.3\n"
            "C\t0\t0.8\n"
            "D\t1\t3.1\n"
        )
        tree_tips = ["A", "B", "C", "D"]
        t1, t2, names = ThresholdModel._parse_multi_trait_file(
            str(trait_file), "t1", "t2", "discrete", "continuous", tree_tips
        )
        assert set(names) == {"A", "B", "C", "D"}
        assert t1["A"] == 0.0
        assert t1["B"] == 1.0
        assert abs(t2["A"] - 1.5) < 1e-10

    def test_comments_blanks_and_extra_columns(self, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "taxon\tt1\tt2\n"
            "# comment after header\n"
            "\n"
            "A\t0\t1.5\textra\tignored\n"
            "B\t1\t2.3\textra\tignored\n"
            "C\t0\t0.8\textra\tignored\n"
            "D\t1\t3.1\textra\tignored\n"
        )
        tree_tips = ["A", "B", "C", "D"]

        t1, t2, names = ThresholdModel._parse_multi_trait_file(
            str(trait_file), "t1", "t2", "discrete", "continuous", tree_tips
        )

        assert names == ["A", "B", "C", "D"]
        assert t1 == {"A": 0.0, "B": 1.0, "C": 0.0, "D": 1.0}
        assert t2["D"] == 3.1

    def test_all_shared_taxa_emit_no_warnings(self, tmp_path, capsys):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "taxon\tt1\tt2\tnote\n"
            "A\t0\t1.5\tx\n"
            "B\t1\t2.3\tx\n"
            "C\t0\t0.8\tx\n"
            "D\t1\t3.1\tx\n"
        )

        t1, t2, names = ThresholdModel._parse_multi_trait_file(
            str(trait_file), "t1", "t2", "discrete", "continuous", ["D", "B", "A", "C"]
        )

        assert names == ["A", "B", "C", "D"]
        assert t1 == {"A": 0.0, "B": 1.0, "C": 0.0, "D": 1.0}
        assert t2 == {"A": 1.5, "B": 2.3, "C": 0.8, "D": 3.1}
        assert capsys.readouterr().err == ""

    def test_first_two_trait_columns_can_be_requested_in_reverse_order(self, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "taxon\tt1\tt2\tnote\n"
            "A\t0\t1.5\tx\n"
            "B\t1\t2.3\tx\n"
            "C\t0\t0.8\tx\n"
        )

        t2, t1, names = ThresholdModel._parse_multi_trait_file(
            str(trait_file), "t2", "t1", "continuous", "discrete", ["A", "B", "C"]
        )

        assert names == ["A", "B", "C"]
        assert t2 == {"A": 1.5, "B": 2.3, "C": 0.8}
        assert t1 == {"A": 0.0, "B": 1.0, "C": 0.0}

    def test_short_row_raises(self, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "taxon\tt1\tt2\n"
            "A\t0\t1.5\n"
            "B\t1\n"
            "C\t0\t0.8\n"
        )

        with pytest.raises(SystemExit):
            ThresholdModel._parse_multi_trait_file(
                str(trait_file), "t1", "t2",
                "discrete", "continuous", ["A", "B", "C"]
            )

    def test_missing_file_raises(self):
        with pytest.raises(SystemExit):
            ThresholdModel._parse_multi_trait_file(
                "/nonexistent/file.tsv", "t1", "t2",
                "discrete", "continuous", ["A"]
            )

    def test_column_not_found_raises(self, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text("taxon\tt1\tt2\nA\t0\t1.0\n")
        with pytest.raises(SystemExit):
            ThresholdModel._parse_multi_trait_file(
                str(trait_file), "t1", "missing_col",
                "discrete", "continuous", ["A"]
            )

    def test_non_binary_raises(self, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "taxon\tt1\tt2\n"
            "A\t0\t1.0\n"
            "B\t2\t2.0\n"
            "C\t1\t3.0\n"
        )
        with pytest.raises(SystemExit):
            ThresholdModel._parse_multi_trait_file(
                str(trait_file), "t1", "t2",
                "discrete", "continuous", ["A", "B", "C"]
            )

    def test_taxon_mismatch_warning(self, tmp_path, capsys):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "taxon\tt1\tt2\n"
            "A\t0\t1.0\n"
            "B\t1\t2.0\n"
            "C\t0\t3.0\n"
            "D\t1\t4.0\n"
            "E\t0\t5.0\n"
        )
        tree_tips = ["A", "B", "C", "F"]
        t1, t2, names = ThresholdModel._parse_multi_trait_file(
            str(trait_file), "t1", "t2",
            "discrete", "continuous", tree_tips
        )
        err = capsys.readouterr().err
        assert "Warning" in err
        assert set(names) == {"A", "B", "C"}

    def test_too_few_shared_taxa_raises(self, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text("taxon\tt1\tt2\nA\t0\t1.0\nB\t1\t2.0\n")
        with pytest.raises(SystemExit):
            ThresholdModel._parse_multi_trait_file(
                str(trait_file), "t1", "t2",
                "discrete", "continuous", ["A", "B"]
            )


class TestBuildVCVMatrix:
    def test_symmetric(self):
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        vcv = ThresholdModel._build_vcv_matrix(tree, names)
        np.testing.assert_array_almost_equal(vcv, vcv.T)

    def test_correct_shape(self):
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        vcv = ThresholdModel._build_vcv_matrix(tree, names)
        assert vcv.shape == (4, 4)

    def test_diagonal_root_to_tip(self):
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        vcv = ThresholdModel._build_vcv_matrix(tree, names)
        for i, name in enumerate(names):
            expected = tree.distance(tree.root, name)
            assert abs(vcv[i, i] - expected) < 1e-10

    def test_prune_tree_to_taxa_uses_direct_batch_pruning(self, monkeypatch):
        tree = _make_tree()
        original_get_terminals = TreeMixin.get_terminals

        def fail_traversal(*_args, **_kwargs):
            raise AssertionError("generic terminal traversal should not be used")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)
        monkeypatch.setattr(TreeMixin, "prune", fail_traversal)

        ThresholdModel._prune_tree_to_taxa(tree, {"A", "C"})

        assert {tip.name for tip in original_get_terminals(tree)} == {"A", "C"}


class TestInitializeLiabilities:
    def test_discrete_signs_correct(self):
        rng = np.random.default_rng(42)
        trait_values = {"A": 0.0, "B": 1.0, "C": 0.0, "D": 1.0}
        names = ["A", "B", "C", "D"]
        liabs = ThresholdModel._initialize_liabilities(
            trait_values, "discrete", names, rng
        )
        assert liabs[0] < 0  # A: state 0
        assert liabs[1] > 0  # B: state 1
        assert liabs[2] < 0  # C: state 0
        assert liabs[3] > 0  # D: state 1

    def test_discrete_initialization_matches_scalar_truncnorm_stream(self):
        from scipy.stats import truncnorm

        trait_values = {
            f"t{i}": float(state)
            for i, state in enumerate([0, 1, 0, 1, 1, 0, 0, 1])
        }
        names = list(trait_values)

        expected = []
        rng_ref = np.random.default_rng(123)
        for name in names:
            if int(trait_values[name]) == 0:
                expected.append(
                    truncnorm.rvs(
                        -np.inf,
                        0,
                        loc=0,
                        scale=1,
                        random_state=rng_ref,
                    )
                )
            else:
                expected.append(
                    truncnorm.rvs(
                        0,
                        np.inf,
                        loc=0,
                        scale=1,
                        random_state=rng_ref,
                    )
                )

        observed = ThresholdModel._initialize_liabilities(
            trait_values,
            "discrete",
            names,
            np.random.default_rng(123),
        )

        np.testing.assert_allclose(observed, expected)

    def test_continuous_passthrough(self):
        rng = np.random.default_rng(42)
        trait_values = {"A": 1.5, "B": -0.3, "C": 2.1, "D": 0.7}
        names = ["A", "B", "C", "D"]
        liabs = ThresholdModel._initialize_liabilities(
            trait_values, "continuous", names, rng
        )
        assert abs(liabs[0] - 1.5) < 1e-10
        assert abs(liabs[1] - (-0.3)) < 1e-10


class TestLogLikelihoodBivariateBM:
    def test_sufficient_stats_match_explicit_residual_quadratics(self):
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        C = ThresholdModel._build_vcv_matrix(tree, names)
        C_inv = np.linalg.inv(C)
        _, logdet_C = np.linalg.slogdet(C)
        n = len(names)

        x1 = np.array([1.0, 1.5, 0.5, 2.0])
        x2 = np.array([0.3, -0.5, 1.0, 0.8])
        sigma2_1 = 1.3
        sigma2_2 = 0.8
        r = 0.35
        a1 = 0.9
        a2 = 0.2

        stats = ThresholdModel._bivariate_quadratic_stats(x1, x2, C_inv)
        cached_ll = ThresholdModel._log_likelihood_bivariate_bm_from_stats(
            stats, sigma2_1, sigma2_2, r, a1, a2, logdet_C, n
        )

        cov12 = r * np.sqrt(sigma2_1 * sigma2_2)
        Sigma = np.array([[sigma2_1, cov12], [cov12, sigma2_2]])
        Sigma_inv = np.linalg.inv(Sigma)
        residuals = np.column_stack([x1 - a1, x2 - a2])
        quad = 0.0
        for i in range(n):
            for j in range(n):
                quad += C_inv[i, j] * float(
                    residuals[i] @ Sigma_inv @ residuals[j]
                )
        _, logdet_Sigma = np.linalg.slogdet(Sigma)
        explicit_ll = -0.5 * (
            2 * n * np.log(2 * np.pi)
            + n * logdet_Sigma
            + 2.0 * logdet_C
            + quad
        )

        assert cached_ll == pytest.approx(explicit_ll)

    def test_sufficient_stats_from_products_match_full_stats(self, monkeypatch):
        C = np.array([
            [1.0, 0.2, 0.1],
            [0.2, 1.0, 0.3],
            [0.1, 0.3, 1.0],
        ])
        C_inv = np.linalg.inv(C)
        x1 = np.array([0.4, -0.2, 1.1])
        x2 = np.array([1.5, 0.7, -0.3])

        expected = ThresholdModel._bivariate_quadratic_stats(x1, x2, C_inv)
        one_C_one = float(np.sum(C_inv))
        monkeypatch.setattr(
            threshold_model_module.np,
            "sum",
            lambda *args, **kwargs: pytest.fail(
                "cached bivariate stats should use vector.sum"
            ),
        )
        observed = ThresholdModel._bivariate_quadratic_stats_from_products(
            x1,
            x2,
            C_inv @ x1,
            C_inv @ x2,
            one_C_one,
        )

        np.testing.assert_allclose(observed, expected)

    def test_r_zero_equals_sum_univariate(self):
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        C = ThresholdModel._build_vcv_matrix(tree, names)
        C_inv = np.linalg.inv(C)
        _, logdet_C = np.linalg.slogdet(C)
        n = len(names)

        x1 = np.array([1.0, 1.5, 0.5, 2.0])
        x2 = np.array([0.3, -0.5, 1.0, 0.8])

        ll_biv = ThresholdModel._log_likelihood_bivariate_bm(
            x1, x2, C_inv, 1.0, 1.0, 0.0, 1.0, 0.5, logdet_C, n
        )

        # At r=0, bivariate likelihood should decompose to two univariates
        r1 = x1 - 1.0
        r2 = x2 - 0.5
        ll1 = -0.5 * (n * np.log(2 * np.pi) + logdet_C + n * np.log(1.0) + float(r1 @ C_inv @ r1))
        ll2 = -0.5 * (n * np.log(2 * np.pi) + logdet_C + n * np.log(1.0) + float(r2 @ C_inv @ r2))

        assert abs(ll_biv - (ll1 + ll2)) < 1e-6

    def test_positive_r_data_higher_ll(self):
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        C = ThresholdModel._build_vcv_matrix(tree, names)
        C_inv = np.linalg.inv(C)
        _, logdet_C = np.linalg.slogdet(C)
        n = len(names)

        # Positively correlated data
        x1 = np.array([1.0, 2.0, 3.0, 4.0])
        x2 = np.array([1.1, 2.1, 2.9, 3.8])

        ll_pos = ThresholdModel._log_likelihood_bivariate_bm(
            x1, x2, C_inv, 1.0, 1.0, 0.8, 2.5, 2.5, logdet_C, n
        )
        ll_neg = ThresholdModel._log_likelihood_bivariate_bm(
            x1, x2, C_inv, 1.0, 1.0, -0.8, 2.5, 2.5, logdet_C, n
        )

        assert ll_pos > ll_neg


class TestSampleLiabilitiesGibbs:
    def test_inverse_cdf_truncated_normal_respects_bounds(self):
        rng = np.random.default_rng(42)
        lower_samples = np.array([
            ThresholdModel._sample_truncated_normal(0.25, 1.2, -np.inf, 0.0, rng)
            for _ in range(100)
        ])
        upper_samples = np.array([
            ThresholdModel._sample_truncated_normal(0.25, 1.2, 0.0, np.inf, rng)
            for _ in range(100)
        ])

        assert np.all(lower_samples <= 0.0)
        assert np.all(upper_samples >= 0.0)

    def test_one_sided_truncated_normal_matches_inverse_cdf(self):
        def reference(mu, sd, lower, upper, rng):
            from scipy.special import ndtr, ndtri

            lower_z = (lower - mu) / sd
            upper_z = (upper - mu) / sd
            lower_cdf = 0.0 if np.isneginf(lower_z) else ndtr(lower_z)
            upper_cdf = 1.0 if np.isposinf(upper_z) else ndtr(upper_z)
            u = lower_cdf + rng.random() * (upper_cdf - lower_cdf)
            u = np.clip(u, np.finfo(float).tiny, 1.0 - np.finfo(float).eps)
            return float(mu + sd * ndtri(u))

        cases = [
            (0.25, 1.2, -np.inf, 0.0),
            (-0.5, 0.7, 0.0, np.inf),
        ]
        for mu, sd, lower, upper in cases:
            rng_fast = np.random.default_rng(123)
            rng_ref = np.random.default_rng(123)
            fast = [
                ThresholdModel._sample_truncated_normal(
                    mu, sd, lower, upper, rng_fast
                )
                for _ in range(50)
            ]
            expected = [
                reference(mu, sd, lower, upper, rng_ref)
                for _ in range(50)
            ]

            np.testing.assert_allclose(fast, expected)

    def test_truncated_normal_uses_cached_float_bounds(self, monkeypatch):
        def fail_finfo(*_args, **_kwargs):
            raise AssertionError("np.finfo should not be called per draw")

        monkeypatch.setattr(
            "phykit.services.tree.threshold_model.np.finfo",
            fail_finfo,
        )

        sample = ThresholdModel._sample_truncated_normal(
            0.25,
            1.2,
            -np.inf,
            0.0,
            np.random.default_rng(123),
        )

        assert sample <= 0.0

    def test_inverse_cdf_truncated_normal_does_not_import_scipy_stats(self, monkeypatch):
        original_import = __import__

        def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name == "scipy.stats" or name.startswith("scipy.stats."):
                raise AssertionError("inverse-CDF truncated normal should not import scipy.stats")
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr("builtins.__import__", fake_import)

        sample = ThresholdModel._sample_truncated_normal(
            0.25,
            1.2,
            -np.inf,
            0.0,
            np.random.default_rng(123),
        )

        assert sample <= 0.0

    def test_one_sided_truncated_normal_does_not_use_cached_ndtr(self, monkeypatch):
        def fail_ndtr(*_args, **_kwargs):
            raise AssertionError("one-sided truncated normal should use math CDF")

        monkeypatch.setattr(threshold_model_module, "_NDTR", fail_ndtr)

        sample = ThresholdModel._sample_truncated_normal(
            0.25,
            1.2,
            -np.inf,
            0.0,
            np.random.default_rng(123),
        )

        assert sample <= 0.0

    def test_inverse_cdf_helpers_are_cached_after_first_draw(self, monkeypatch):
        threshold_model_module._NDTR = None
        threshold_model_module._NDTRI = None
        original_import = builtins.__import__
        special_imports = 0

        def tracking_import(name, globals=None, locals=None, fromlist=(), level=0):
            nonlocal special_imports
            if name == "scipy.special":
                special_imports += 1
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", tracking_import)
        rng = np.random.default_rng(123)

        for _ in range(5):
            ThresholdModel._sample_truncated_normal(
                0.25,
                1.2,
                -np.inf,
                0.0,
                rng,
            )

        assert special_imports == 1
        assert threshold_model_module._NDTR is None
        assert threshold_model_module._NDTRI is not None

    def test_liabilities_respect_threshold(self):
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        C = ThresholdModel._build_vcv_matrix(tree, names)
        C_inv = np.linalg.inv(C)
        rng = np.random.default_rng(42)

        states = np.array([0, 1, 0, 1])
        liabilities = np.array([-0.5, 0.5, -0.3, 0.8])
        other_liabs = np.array([1.0, 2.0, 0.5, 1.5])

        new_liabs = ThresholdModel._sample_liabilities_gibbs(
            liabilities, states, C, C_inv, 1.0, 0.0,
            other_liabs, 1.0, 0.0, 1.0, rng
        )

        # State 0 must be negative, state 1 must be positive
        assert new_liabs[0] < 0
        assert new_liabs[1] > 0
        assert new_liabs[2] < 0
        assert new_liabs[3] > 0

    def test_deterministic_with_seed(self):
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        C = ThresholdModel._build_vcv_matrix(tree, names)
        C_inv = np.linalg.inv(C)

        states = np.array([0, 1, 0, 1])
        liabilities = np.array([-0.5, 0.5, -0.3, 0.8])
        other_liabs = np.array([1.0, 2.0, 0.5, 1.5])

        rng1 = np.random.default_rng(123)
        result1 = ThresholdModel._sample_liabilities_gibbs(
            liabilities, states, C, C_inv, 1.0, 0.0,
            other_liabs, 1.0, 0.0, 1.0, rng1
        )

        rng2 = np.random.default_rng(123)
        result2 = ThresholdModel._sample_liabilities_gibbs(
            liabilities, states, C, C_inv, 1.0, 0.0,
            other_liabs, 1.0, 0.0, 1.0, rng2
        )

        np.testing.assert_array_almost_equal(result1, result2)

    def test_gibbs_sampler_matches_materialized_residual_reference(self):
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        C = ThresholdModel._build_vcv_matrix(tree, names)
        C_inv = np.linalg.inv(C)

        states = np.array([0, 1, 0, 1])
        liabilities = np.array([-0.5, 0.5, -0.3, 0.8])
        other_liabs = np.array([1.0, 2.0, 0.5, 1.5])
        sigma2 = 1.2
        a = 0.1

        def reference(rng):
            new_liabilities = liabilities.copy()
            for i in range(len(liabilities)):
                c_ii_inv = C_inv[i, i]
                var_cond = sigma2 / c_ii_inv
                residuals = new_liabilities - a
                row_sum = float(C_inv[i, :] @ residuals) - c_ii_inv * residuals[i]
                mu_cond = a - row_sum / c_ii_inv
                sd_cond = np.sqrt(var_cond)
                if states[i] == 0:
                    new_liabilities[i] = ThresholdModel._sample_truncated_normal(
                        mu_cond, sd_cond, -np.inf, 0.0, rng
                    )
                else:
                    new_liabilities[i] = ThresholdModel._sample_truncated_normal(
                        mu_cond, sd_cond, 0.0, np.inf, rng
                    )
            return new_liabilities

        observed = ThresholdModel._sample_liabilities_gibbs(
            liabilities,
            states,
            C,
            C_inv,
            sigma2,
            0.0,
            other_liabs,
            1.0,
            a,
            1.0,
            np.random.default_rng(123),
        )
        expected = reference(np.random.default_rng(123))

        np.testing.assert_allclose(observed, expected)

    def test_gibbs_context_matches_default_sampler(self):
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        C = ThresholdModel._build_vcv_matrix(tree, names)
        C_inv = np.linalg.inv(C)

        states = np.array([0, 1, 0, 1])
        liabilities = np.array([-0.5, 0.5, -0.3, 0.8])
        other_liabs = np.array([1.0, 2.0, 0.5, 1.5])
        sigma2 = 1.2
        a = 0.1
        context = ThresholdModel._prepare_gibbs_context(C_inv, sigma2)

        expected = ThresholdModel._sample_liabilities_gibbs(
            liabilities,
            states,
            C,
            C_inv,
            sigma2,
            0.0,
            other_liabs,
            1.0,
            a,
            1.0,
            np.random.default_rng(123),
        )
        observed = ThresholdModel._sample_liabilities_gibbs(
            liabilities,
            states,
            C,
            C_inv,
            sigma2,
            0.0,
            other_liabs,
            1.0,
            a,
            1.0,
            np.random.default_rng(123),
            context,
        )

        np.testing.assert_array_equal(observed, expected)

    def test_gibbs_context_skips_per_tip_numpy_finite_checks(self, monkeypatch):
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        C = ThresholdModel._build_vcv_matrix(tree, names)
        C_inv = np.linalg.inv(C)

        states = np.array([0, 1, 0, 1])
        liabilities = np.array([-0.5, 0.5, -0.3, 0.8])
        other_liabs = np.array([1.0, 2.0, 0.5, 1.5])
        context = ThresholdModel._prepare_gibbs_context(C_inv, 1.2)

        def fail_isfinite(*_args, **_kwargs):
            raise AssertionError("context path should not check finite values per tip")

        monkeypatch.setattr(
            "phykit.services.tree.threshold_model.np.isfinite",
            fail_isfinite,
        )

        observed = ThresholdModel._sample_liabilities_gibbs(
            liabilities,
            states,
            C,
            C_inv,
            1.2,
            0.0,
            other_liabs,
            1.0,
            0.1,
            1.0,
            np.random.default_rng(123),
            context,
        )

        assert observed.shape == liabilities.shape

    def test_gibbs_sampler_uses_one_sided_inverse_cdf_path(self, monkeypatch):
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        C = ThresholdModel._build_vcv_matrix(tree, names)
        C_inv = np.linalg.inv(C)

        states = np.array([0, 1, 0, 1])
        liabilities = np.array([-0.5, 0.5, -0.3, 0.8])
        other_liabs = np.array([1.0, 2.0, 0.5, 1.5])

        def fail_generic_sampler(*_args, **_kwargs):
            raise AssertionError("Gibbs should avoid the generic sampler dispatcher")

        monkeypatch.setattr(
            ThresholdModel,
            "_sample_truncated_normal",
            staticmethod(fail_generic_sampler),
        )

        observed = ThresholdModel._sample_liabilities_gibbs(
            liabilities,
            states,
            C,
            C_inv,
            1.2,
            0.0,
            other_liabs,
            1.0,
            0.1,
            1.0,
            np.random.default_rng(123),
        )

        assert np.all(observed[states == 0] <= 0.0)
        assert np.all(observed[states == 1] >= 0.0)


class TestComputeHPD:
    def test_known_normal_samples(self):
        rng = np.random.default_rng(42)
        samples = rng.normal(0, 1, 10000)
        lo, hi = ThresholdModel._compute_hpd(samples, 0.95)
        # 95% HPD of standard normal should be approximately [-1.96, 1.96]
        assert lo < -1.5
        assert hi > 1.5
        assert (hi - lo) < 5.0  # Not absurdly wide

    def test_95_coverage(self):
        rng = np.random.default_rng(42)
        samples = rng.normal(5, 2, 10000)
        lo, hi = ThresholdModel._compute_hpd(samples, 0.95)
        in_interval = np.sum((samples >= lo) & (samples <= hi))
        coverage = in_interval / len(samples)
        assert coverage >= 0.94

    def test_hpd_from_sorted_matches_public_helper(self):
        rng = np.random.default_rng(20260628)
        samples = rng.normal(size=101)
        sorted_samples = np.sort(samples)

        assert ThresholdModel._compute_hpd_from_sorted(sorted_samples) == pytest.approx(
            ThresholdModel._compute_hpd(samples)
        )

    def test_summarize_posterior_reuses_sorted_samples_for_median(self, monkeypatch):
        rng = np.random.default_rng(20260628)
        samples = rng.normal(size=101)
        mcmc_result = {
            "r": samples,
            "sigma2_1": samples + 1.0,
            "sigma2_2": samples + 2.0,
            "a1": samples - 1.0,
            "a2": samples - 2.0,
            "acceptance_rates": {"r": 0.25},
        }

        def fail_median(*_args, **_kwargs):
            raise AssertionError("summary should reuse the HPD sorted samples")

        monkeypatch.setattr(threshold_model_module.np, "median", fail_median)

        summary = ThresholdModel._summarize_posterior(mcmc_result)

        assert summary["r"]["mean"] == pytest.approx(float(np.mean(samples)))
        assert summary["r"]["median"] == pytest.approx(float(np.median(samples)))
        assert summary["acceptance_rates"] == {"r": 0.25}

    def test_summarize_posterior_small_traces_use_array_mean(self, monkeypatch):
        rng = np.random.default_rng(20260701)
        samples = rng.normal(size=101)
        mcmc_result = {
            "r": samples,
            "sigma2_1": samples + 1.0,
            "sigma2_2": samples + 2.0,
            "a1": samples - 1.0,
            "a2": samples - 2.0,
            "acceptance_rates": {"r": 0.25},
        }

        def fail_mean(*_args, **_kwargs):
            raise AssertionError("small posterior traces should use ndarray.mean")

        monkeypatch.setattr(threshold_model_module.np, "mean", fail_mean)

        summary = ThresholdModel._summarize_posterior(mcmc_result)

        assert summary["r"]["mean"] == pytest.approx(float(samples.mean()))


class TestRunMCMC:
    def test_run_mcmc_mean_diagonal_uses_trace(self, monkeypatch):
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        C = ThresholdModel._build_vcv_matrix(tree, names)
        rng = np.random.default_rng(42)
        t1 = {"A": 0.2, "B": 1.0, "C": 0.4, "D": 1.2}
        t2 = {"A": 1.0, "B": 2.0, "C": 0.5, "D": 1.8}

        def fail_diag(*_args, **_kwargs):
            raise AssertionError("MCMC setup should use np.trace instead of np.diag")

        def fail_reduction(*_args, **_kwargs):
            raise AssertionError("continuous MCMC setup should use ndarray reductions")

        monkeypatch.setattr(
            ThresholdModel,
            "_vcv_inverse_and_logdet",
            staticmethod(lambda matrix: (np.linalg.inv(matrix), 0.0)),
        )
        monkeypatch.setattr(threshold_model_module.np, "diag", fail_diag)
        monkeypatch.setattr(threshold_model_module.np, "var", fail_reduction)
        monkeypatch.setattr(threshold_model_module.np, "mean", fail_reduction)

        result = ThresholdModel._run_mcmc(
            t1, t2, "continuous", "continuous",
            names, C, 30, 10, 0.2, rng,
        )

        assert "sigma2_2" in result

    def test_expected_keys(self):
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        C = ThresholdModel._build_vcv_matrix(tree, names)
        rng = np.random.default_rng(42)

        t1 = {"A": 0.0, "B": 1.0, "C": 0.0, "D": 1.0}
        t2 = {"A": 1.0, "B": 2.0, "C": 0.5, "D": 1.8}

        result = ThresholdModel._run_mcmc(
            t1, t2, "discrete", "continuous",
            names, C, 1000, 10, 0.2, rng
        )

        assert "r" in result
        assert "sigma2_1" in result
        assert "sigma2_2" in result
        assert "a1" in result
        assert "a2" in result
        assert "acceptance_rates" in result

    def test_sample_count(self):
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        C = ThresholdModel._build_vcv_matrix(tree, names)
        rng = np.random.default_rng(42)

        t1 = {"A": 0.0, "B": 1.0, "C": 0.0, "D": 1.0}
        t2 = {"A": 1.0, "B": 2.0, "C": 0.5, "D": 1.8}

        ngen = 2000
        sample = 10
        burnin = 0.5
        result = ThresholdModel._run_mcmc(
            t1, t2, "discrete", "continuous",
            names, C, ngen, sample, burnin, rng
        )

        expected = (ngen - int(ngen * burnin)) // sample
        assert len(result["r"]) == expected

    def test_r_in_bounds(self):
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        C = ThresholdModel._build_vcv_matrix(tree, names)
        rng = np.random.default_rng(42)

        t1 = {"A": 0.0, "B": 1.0, "C": 0.0, "D": 1.0}
        t2 = {"A": 1.0, "B": 2.0, "C": 0.5, "D": 1.8}

        result = ThresholdModel._run_mcmc(
            t1, t2, "discrete", "continuous",
            names, C, 2000, 10, 0.2, rng
        )

        assert np.all(result["r"] >= -1.0)
        assert np.all(result["r"] <= 1.0)

    def test_sigma2_positive(self):
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        C = ThresholdModel._build_vcv_matrix(tree, names)
        rng = np.random.default_rng(42)

        t1 = {"A": 0.0, "B": 1.0, "C": 0.0, "D": 1.0}
        t2 = {"A": 1.0, "B": 2.0, "C": 0.5, "D": 1.8}

        result = ThresholdModel._run_mcmc(
            t1, t2, "discrete", "continuous",
            names, C, 2000, 10, 0.2, rng
        )

        assert np.all(result["sigma2_1"] > 0)
        assert np.all(result["sigma2_2"] > 0)

    def test_reproducible_with_seed(self):
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        C = ThresholdModel._build_vcv_matrix(tree, names)

        t1 = {"A": 0.0, "B": 1.0, "C": 0.0, "D": 1.0}
        t2 = {"A": 1.0, "B": 2.0, "C": 0.5, "D": 1.8}

        rng1 = np.random.default_rng(42)
        result1 = ThresholdModel._run_mcmc(
            t1, t2, "discrete", "continuous",
            names, C, 1000, 10, 0.2, rng1
        )

        rng2 = np.random.default_rng(42)
        result2 = ThresholdModel._run_mcmc(
            t1, t2, "discrete", "continuous",
            names, C, 1000, 10, 0.2, rng2
        )

        np.testing.assert_array_almost_equal(result1["r"], result2["r"])

    def test_adaptive_tuning_keeps_cached_proposal_scales_valid(self):
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        C = ThresholdModel._build_vcv_matrix(tree, names)
        rng = np.random.default_rng(42)

        t1 = {"A": 1.0, "B": 2.0, "C": 0.5, "D": 1.8}
        t2 = {"A": 1.1, "B": 2.1, "C": 0.6, "D": 1.9}

        result = ThresholdModel._run_mcmc(
            t1, t2, "continuous", "continuous",
            names, C, 300, 10, 0.8, rng
        )

        for param in ("r", "sigma2_1", "sigma2_2", "a1", "a2"):
            assert np.all(np.isfinite(result[param]))
        assert set(result["acceptance_rates"]) == {
            "log_sigma2_1",
            "log_sigma2_2",
            "r",
            "a1",
            "a2",
        }

    def test_continuous_initialization_uses_array_reductions(self, monkeypatch):
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        C = ThresholdModel._build_vcv_matrix(tree, names)
        rng = np.random.default_rng(42)

        t1 = {"A": 1.0, "B": 2.0, "C": 0.5, "D": 1.8}
        t2 = {"A": 1.1, "B": 2.1, "C": 0.6, "D": 1.9}

        def fail_reduction(*_args, **_kwargs):
            raise AssertionError("continuous initialization should use ndarray reductions")

        monkeypatch.setattr(threshold_model_module.np, "var", fail_reduction)
        monkeypatch.setattr(threshold_model_module.np, "mean", fail_reduction)

        result = ThresholdModel._run_mcmc(
            t1, t2, "continuous", "continuous",
            names, C, 2, 1, 0.0, rng
        )

        assert len(result["sigma2_1"]) == 2
        assert len(result["a1"]) == 2

    def test_continuous_continuous(self):
        """Both continuous: should skip Gibbs, just run MH."""
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        C = ThresholdModel._build_vcv_matrix(tree, names)
        rng = np.random.default_rng(42)

        t1 = {"A": 1.0, "B": 2.0, "C": 0.5, "D": 1.8}
        t2 = {"A": 1.1, "B": 2.1, "C": 0.6, "D": 1.9}

        result = ThresholdModel._run_mcmc(
            t1, t2, "continuous", "continuous",
            names, C, 1000, 10, 0.2, rng
        )

        assert len(result["r"]) > 0
        assert np.all(result["sigma2_1"] > 0)


class TestRun:
    @staticmethod
    def _make_args():
        return Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            traits="habitat,body_mass",
            types="discrete,continuous",
            ngen=0,
            sample=1,
            burnin=0.0,
            seed=42,
            plot_output=None,
            json=False,
        )

    @staticmethod
    def _stub_run_tail(monkeypatch, svc, ordered_names, captured):
        trait1 = {name: float(i & 1) for i, name in enumerate(ordered_names)}
        trait2 = {name: float(i) for i, name in enumerate(ordered_names)}

        monkeypatch.setattr(
            ThresholdModel,
            "_parse_multi_trait_file",
            staticmethod(
                lambda *_args, **_kwargs: (trait1, trait2, list(ordered_names))
            ),
        )

        def build_vcv(tree, names):
            captured["tree"] = tree
            captured["names"] = list(names)
            return np.eye(len(names))

        monkeypatch.setattr(svc, "_build_vcv_matrix", build_vcv)
        monkeypatch.setattr(svc, "_run_mcmc", lambda *_args, **_kwargs: {})
        monkeypatch.setattr(svc, "_summarize_posterior", lambda *_args, **_kwargs: {})
        monkeypatch.setattr(svc, "_output_text", lambda *_args, **_kwargs: None)

    def test_run_all_tips_present_uses_read_only_tree_without_copy_or_prune(
        self, monkeypatch
    ):
        tree = _make_tree()
        svc = ThresholdModel(self._make_args())
        captured = {}
        ordered_names = ["A", "B", "C", "D"]

        monkeypatch.setattr(svc, "read_tree_file_unmodified", lambda: tree)
        monkeypatch.setattr(
            svc,
            "read_tree_file",
            lambda: (_ for _ in ()).throw(
                AssertionError("mutable tree read should not be used")
            ),
        )
        monkeypatch.setattr(
            svc,
            "_fast_copy",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("tree should not be copied")
            ),
        )
        monkeypatch.setattr(
            ThresholdModel,
            "_prune_tree_to_taxa",
            staticmethod(
                lambda *_args, **_kwargs: (_ for _ in ()).throw(
                    AssertionError("all-shared data should not prune")
                )
            ),
        )
        self._stub_run_tail(monkeypatch, svc, ordered_names, captured)

        svc.run()

        assert captured["tree"] is tree
        assert captured["names"] == ordered_names
        assert {tip.name for tip in tree.get_terminals()} == {"A", "B", "C", "D"}

    def test_run_missing_trait_taxa_copies_before_pruning(self, monkeypatch):
        tree = _make_tree()
        svc = ThresholdModel(self._make_args())
        original_fast_copy = svc._fast_copy
        copied_trees = []
        captured = {}
        ordered_names = ["A", "B", "C"]

        def copy_spy(tree_to_copy):
            copied_tree = original_fast_copy(tree_to_copy)
            copied_trees.append(copied_tree)
            return copied_tree

        monkeypatch.setattr(svc, "read_tree_file_unmodified", lambda: tree)
        monkeypatch.setattr(svc, "_fast_copy", copy_spy)
        self._stub_run_tail(monkeypatch, svc, ordered_names, captured)

        svc.run()

        assert len(copied_trees) == 1
        assert captured["tree"] is copied_trees[0]
        assert {tip.name for tip in tree.get_terminals()} == {"A", "B", "C", "D"}
        assert {tip.name for tip in copied_trees[0].get_terminals()} == {
            "A",
            "B",
            "C",
        }

    def test_text_output(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            traits="habitat,body_mass",
            types="discrete,continuous",
            ngen=2000,
            sample=10,
            burnin=0.2,
            seed=42,
            plot_output=None,
            json=False,
        )
        svc = ThresholdModel(args)
        svc.run()
        out, _ = capsys.readouterr()
        assert "Trait 1: habitat" in out
        assert "Trait 2: body_mass" in out
        assert "Posterior correlation (r):" in out
        assert "Posterior sigma2_1:" in out
        assert "Acceptance rates:" in out

    def test_output_text_batches_summary(self, mocker):
        svc = ThresholdModel(self._make_args())
        svc.ngen = 2000
        svc.sample = 10
        svc.burnin = 0.2
        summary = {
            "r": {"mean": 0.12345, "hpd_lower": -0.23456, "hpd_upper": 0.45678},
            "sigma2_1": {
                "mean": 1.23456,
                "hpd_lower": 0.98765,
                "hpd_upper": 1.54321,
            },
            "sigma2_2": {
                "mean": 2.34567,
                "hpd_lower": 1.87654,
                "hpd_upper": 2.76543,
            },
            "acceptance_rates": {
                "r": 0.1234,
                "log_sigma2_1": 0.2345,
                "log_sigma2_2": 0.3456,
                "a1": 0.4567,
                "a2": 0.5678,
            },
        }
        printed = mocker.patch("builtins.print")

        svc._output_text(summary)

        printed.assert_called_once_with(
            "Trait 1: habitat (discrete, 2 states: 0, 1)\n"
            "Trait 2: body_mass (continuous)\n"
            "MCMC: 2000 generations, sampled every 10, burn-in 20%\n"
            "---\n"
            "Posterior correlation (r): 0.1235 "
            "(95% HPD: -0.235, 0.457)\n"
            "Posterior sigma2_1: 1.2346 (95% HPD: 0.988, 1.543)\n"
            "Posterior sigma2_2: 2.3457 (95% HPD: 1.877, 2.765)\n"
            "Acceptance rates: r=0.123, sigma2_1=0.234, "
            "sigma2_2=0.346, a1=0.457, a2=0.568"
        )

    def test_json_output(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            traits="habitat,body_mass",
            types="discrete,continuous",
            ngen=2000,
            sample=10,
            burnin=0.2,
            seed=42,
            plot_output=None,
            json=True,
        )
        svc = ThresholdModel(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "metadata" in data
        assert "summary" in data
        assert "posterior_samples" in data
        assert data["metadata"]["trait1"] == "habitat"
        assert "r" in data["summary"]
