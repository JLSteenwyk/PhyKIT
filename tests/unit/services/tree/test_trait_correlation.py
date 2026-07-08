import builtins
import importlib
import json
import os
import subprocess
import sys
import tempfile

import pytest
import numpy as np
from argparse import Namespace
from pathlib import Path

from phykit.services.tree.trait_correlation import TraitCorrelation
import phykit.services.tree.trait_correlation as trait_correlation_module
from phykit.helpers.trait_parsing import parse_multi_trait_file
from phykit.errors import PhykitUserError


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
MULTI_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_multi_traits.tsv")


def test_module_import_does_not_import_scipy_linalg(monkeypatch):
    module_name = "phykit.services.tree.trait_correlation"
    previous = sys.modules.pop(module_name, None)
    original_import = builtins.__import__

    def guarded_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "scipy.linalg" or name.startswith("scipy.linalg."):
            raise AssertionError(
                "trait_correlation module import should not import SciPy linalg"
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


def test_module_import_does_not_import_numpy_or_scipy_linalg():
    code = """
import sys
import phykit.services.tree.trait_correlation as module
assert hasattr(module.np, "__getattr__")
assert callable(module.print_json)
assert callable(module.parse_multi_trait_file)
assert "typing" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
assert "phykit.helpers.trait_parsing" not in sys.modules
assert "numpy" not in sys.modules
assert "scipy.linalg" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_lazy_numpy_caches_resolved_attributes():
    lazy_np = trait_correlation_module._LazyNumpy()

    count_nonzero_attr = lazy_np.count_nonzero

    assert lazy_np.__dict__["count_nonzero"] is count_nonzero_attr
    assert lazy_np.count_nonzero is count_nonzero_attr
    assert lazy_np._module is not None


def test_t_two_tailed_p_values_match_expected_values():
    p_values = trait_correlation_module._t_two_tailed_p_values(
        np.array([0.0, 2.2281388519649385, -2.2281388519649385]),
        10,
    )
    assert p_values == pytest.approx([1.0, 0.05, 0.05])


def test_t_two_tailed_p_values_do_not_import_scipy_stats(monkeypatch):
    original_import = __import__

    def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "scipy.stats" or name.startswith("scipy.stats."):
            raise AssertionError("trait correlation p-values should not import scipy.stats")
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr("builtins.__import__", fake_import)

    p_values = trait_correlation_module._t_two_tailed_p_values(
        np.array([0.0, 2.2281388519649385]),
        10,
    )
    assert p_values == pytest.approx([1.0, 0.05])


def test_t_two_tailed_p_values_cache_scipy_special(monkeypatch):
    trait_correlation_module._STDTR = None
    first = trait_correlation_module._t_two_tailed_p_values(
        np.array([0.0, 2.2281388519649385]),
        10,
    )
    cached_stdtr = trait_correlation_module._STDTR
    assert cached_stdtr is not None

    original_import = __import__

    def fail_special_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "scipy.special" or name.startswith("scipy.special."):
            raise AssertionError("stdtr should be reused after first resolution")
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr("builtins.__import__", fail_special_import)

    second = trait_correlation_module._t_two_tailed_p_values(
        np.array([0.0, 2.2281388519649385]),
        10,
    )
    assert trait_correlation_module._STDTR is cached_stdtr
    assert first == pytest.approx(second)


def test_cholesky_wrappers_cache_scipy_linalg(monkeypatch):
    previous_cho_factor = trait_correlation_module._CHO_FACTOR
    previous_cho_solve = trait_correlation_module._CHO_SOLVE
    trait_correlation_module._CHO_FACTOR = None
    trait_correlation_module._CHO_SOLVE = None

    try:
        factor = trait_correlation_module.cho_factor(
            np.eye(2), lower=True, check_finite=False
        )
        first = trait_correlation_module.cho_solve(
            factor, np.ones(2), check_finite=False
        )
        cached_cho_factor = trait_correlation_module._CHO_FACTOR
        cached_cho_solve = trait_correlation_module._CHO_SOLVE
        assert cached_cho_factor is not None
        assert cached_cho_solve is not None

        original_import = __import__

        def fail_linalg_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name == "scipy.linalg" or name.startswith("scipy.linalg."):
                raise AssertionError(
                    "SciPy linalg functions should be reused after first resolution"
                )
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr("builtins.__import__", fail_linalg_import)

        factor = trait_correlation_module.cho_factor(
            np.eye(2), lower=True, check_finite=False
        )
        second = trait_correlation_module.cho_solve(
            factor, np.ones(2), check_finite=False
        )

        assert trait_correlation_module._CHO_FACTOR is cached_cho_factor
        assert trait_correlation_module._CHO_SOLVE is cached_cho_solve
        assert second == pytest.approx(first)
    finally:
        trait_correlation_module._CHO_FACTOR = previous_cho_factor
        trait_correlation_module._CHO_SOLVE = previous_cho_solve


def _make_args(output_path, **overrides):
    defaults = dict(
        tree=TREE_SIMPLE,
        trait_data=MULTI_TRAITS_FILE,
        output=output_path,
        alpha=0.05,
        cluster=False,
        gene_trees=None,
        json=False,
    )
    defaults.update(overrides)
    return Namespace(**defaults)


def _build_service(output_path, **overrides):
    args = _make_args(output_path, **overrides)
    return TraitCorrelation(args)


class TestTraitCorrelation:
    def test_run_uses_unmodified_tree_read(self, mocker, tmp_path):
        out = str(tmp_path / "corr.png")
        args = _make_args(out)
        svc = TraitCorrelation(args)
        tree = object()
        mocked_read = mocker.patch.object(
            TraitCorrelation,
            "read_tree_file_unmodified",
            return_value=tree,
        )
        mocker.patch.object(
            TraitCorrelation,
            "read_tree_file",
            side_effect=AssertionError("run should not copy cached trees"),
        )
        mocked_validate = mocker.patch.object(TraitCorrelation, "validate_tree")
        mocker.patch.object(
            TraitCorrelation,
            "get_tip_names_from_tree",
            return_value=["a", "b", "c"],
        )
        mocker.patch(
            "phykit.services.tree.trait_correlation.parse_multi_trait_file",
            return_value=(
                ["x", "y", "z"],
                {
                    "a": [1.0, 2.0, 3.0],
                    "b": [2.0, 3.0, 4.0],
                    "c": [3.0, 4.0, 5.0],
                },
            ),
        )
        mocked_vcv = mocker.patch(
            "phykit.services.tree.vcv_utils.build_vcv_matrix",
            return_value=np.eye(3),
        )
        corr = np.eye(3)
        pvals = np.ones((3, 3))
        mocker.patch.object(
            TraitCorrelation,
            "_compute_correlation_matrices",
            return_value=(corr, pvals),
        )
        mocker.patch.object(TraitCorrelation, "_plot_heatmap")
        mocker.patch.object(TraitCorrelation, "_print_text")

        svc.run()

        mocked_read.assert_called_once_with()
        mocked_validate.assert_called_once_with(
            tree,
            min_tips=3,
            require_branch_lengths=True,
            context="trait correlation analysis",
        )
        mocked_vcv.assert_called_once_with(tree, ["a", "b", "c"])

    def test_correlation_matrix_symmetric(self, tmp_path):
        out = str(tmp_path / "corr.png")
        svc = _build_service(out)
        svc.run()

        # Re-run internal logic to inspect the matrix
        corr, pmat, trait_names, n = self._compute_correlation(svc)
        for i in range(len(trait_names)):
            for j in range(len(trait_names)):
                assert abs(corr[i, j] - corr[j, i]) < 1e-10, \
                    f"corr[{i},{j}] != corr[{j},{i}]"

    def test_correlation_diagonal_is_one(self, tmp_path):
        out = str(tmp_path / "corr.png")
        svc = _build_service(out)
        svc.run()

        corr, pmat, trait_names, n = self._compute_correlation(svc)
        for i in range(len(trait_names)):
            assert abs(corr[i, i] - 1.0) < 1e-10

    def test_correlation_range(self, tmp_path):
        out = str(tmp_path / "corr.png")
        svc = _build_service(out)
        svc.run()

        corr, pmat, trait_names, n = self._compute_correlation(svc)
        assert np.all(corr >= -1.0 - 1e-10)
        assert np.all(corr <= 1.0 + 1e-10)

    def test_p_values_diagonal_is_one(self, tmp_path):
        out = str(tmp_path / "corr.png")
        svc = _build_service(out)
        svc.run()

        corr, pmat, trait_names, n = self._compute_correlation(svc)
        for i in range(len(trait_names)):
            assert abs(pmat[i, i] - 1.0) < 1e-10

    def test_p_values_nonnegative(self, tmp_path):
        out = str(tmp_path / "corr.png")
        svc = _build_service(out)
        svc.run()

        corr, pmat, trait_names, n = self._compute_correlation(svc)
        assert np.all(pmat >= 0.0)

    def test_cholesky_correlation_matches_inverse(self, tmp_path):
        rng = np.random.default_rng(20260623)
        out = str(tmp_path / "corr.png")
        svc = _build_service(out)
        n = 16
        p = 4
        A = rng.normal(size=(n, n))
        vcv = A @ A.T + np.eye(n)
        Y = rng.normal(size=(n, p))

        fast_corr, fast_p = svc._compute_correlation_matrices_cholesky(Y, vcv)
        inverse_corr, inverse_p = svc._compute_correlation_matrices_inverse(Y, vcv)

        assert fast_corr == pytest.approx(inverse_corr)
        assert fast_p == pytest.approx(inverse_p)

    def test_cholesky_correlation_uses_single_solve(self, tmp_path, monkeypatch):
        rng = np.random.default_rng(20260628)
        out = str(tmp_path / "corr.png")
        svc = _build_service(out)
        n = 16
        p = 4
        A = rng.normal(size=(n, n))
        vcv = A @ A.T + np.eye(n)
        Y = rng.normal(size=(n, p))
        original_cho_solve = trait_correlation_module.cho_solve
        solve_calls = 0

        def counting_cho_solve(*args, **kwargs):
            nonlocal solve_calls
            solve_calls += 1
            return original_cho_solve(*args, **kwargs)

        monkeypatch.setattr(trait_correlation_module, "cho_solve", counting_cho_solve)

        fast_corr, fast_p = svc._compute_correlation_matrices_cholesky(Y, vcv)
        inverse_corr, inverse_p = svc._compute_correlation_matrices_inverse(Y, vcv)

        assert solve_calls == 1
        assert fast_corr == pytest.approx(inverse_corr)
        assert fast_p == pytest.approx(inverse_p)

    def test_inverse_correlation_matches_scalar_reference(self, tmp_path):
        rng = np.random.default_rng(20260624)
        out = str(tmp_path / "corr.png")
        svc = _build_service(out)
        n = 14
        p = 5
        A = rng.normal(size=(n, n))
        vcv = A @ A.T + np.eye(n)
        Y = rng.normal(size=(n, p))
        observed_corr, observed_p = svc._compute_correlation_matrices_inverse(Y, vcv)

        C_inv = np.linalg.inv(vcv)
        ones = np.ones(n)
        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])
        Y_centered = Y - a_hat
        phylo_cov = (Y_centered.T @ C_inv @ Y_centered) / (n - 1)
        expected_corr, expected_p = svc._correlation_and_p_values(phylo_cov, n, p)

        assert observed_corr == pytest.approx(expected_corr)
        assert observed_p == pytest.approx(expected_p)

    def test_correlation_p_values_match_scalar_reference(
        self, tmp_path, monkeypatch
    ):
        from scipy.stats import t as t_dist

        out = str(tmp_path / "corr.png")
        svc = _build_service(out)
        phylo_cov = np.array(
            [
                [4.0, 2.0, 1.0],
                [2.0, 1.0, 0.5],
                [1.0, 0.5, 9.0],
            ]
        )
        n = 12

        def fail_diag(*_args, **_kwargs):
            raise AssertionError("correlation scaling should use diagonal access")

        monkeypatch.setattr(trait_correlation_module.np, "diag", fail_diag)

        corr, pmat = svc._correlation_and_p_values(phylo_cov, n, 3)

        assert pmat.diagonal().tolist() == [1.0, 1.0, 1.0]
        assert pmat[0, 1] == 0.0
        assert pmat[1, 0] == 0.0
        r = corr[0, 2]
        t_stat = r * np.sqrt((n - 2) / (1 - r ** 2))
        expected_p = 2.0 * t_dist.sf(abs(t_stat), df=n - 2)
        assert pmat[0, 2] == pytest.approx(expected_p)
        assert pmat[2, 0] == pytest.approx(expected_p)

    def test_correlation_p_values_evaluate_upper_triangle_once(
        self, monkeypatch, tmp_path
    ):
        out = str(tmp_path / "corr.png")
        svc = _build_service(out)
        phylo_cov = np.array(
            [
                [4.0, 0.8, 0.4, 0.2],
                [0.8, 9.0, 0.6, 0.3],
                [0.4, 0.6, 16.0, 0.5],
                [0.2, 0.3, 0.5, 25.0],
            ]
        )
        observed_lengths = []

        def fake_p_values(t_stats, df):
            observed_lengths.append(len(t_stats))
            return np.full(len(t_stats), 0.25)

        monkeypatch.setattr(
            trait_correlation_module,
            "_t_two_tailed_p_values",
            fake_p_values,
        )

        _, pmat = svc._correlation_and_p_values(phylo_cov, 12, 4)

        assert observed_lengths == [6]
        assert np.diag(pmat).tolist() == [1.0, 1.0, 1.0, 1.0]
        assert np.allclose(pmat, pmat.T)
        assert pmat[0, 1] == 0.25

    def test_creates_png(self, tmp_path):
        out = str(tmp_path / "corr.png")
        svc = _build_service(out)
        svc.run()
        assert os.path.exists(out)

    def test_cluster_mode(self, tmp_path):
        out = str(tmp_path / "corr_cluster.png")
        svc = _build_service(out, cluster=True)
        svc.run()
        assert os.path.exists(out)

    @pytest.mark.parametrize("cluster", [False, True])
    def test_plot_heatmap_batches_significance_stars(
        self, monkeypatch, tmp_path, cluster
    ):
        pytest.importorskip("matplotlib")
        import matplotlib.axes

        out = str(tmp_path / f"corr_stars_{cluster}.png")
        svc = _build_service(out, cluster=cluster)
        corr = np.eye(4)
        pmat = np.ones((4, 4))
        pmat[0, 1] = pmat[1, 0] = 0.0005
        pmat[0, 2] = pmat[2, 0] = 0.005
        pmat[1, 2] = pmat[2, 1] = 0.02
        trait_names = ["a", "b", "c", "d"]
        original_text = matplotlib.axes.Axes.text
        original_scatter = matplotlib.axes.Axes.scatter
        scatter_sizes = {}

        def fail_star_text(self, *args, **kwargs):
            if len(args) >= 3 and args[2] in {"*", "**", "***"}:
                raise AssertionError("significance stars should be batched")
            return original_text(self, *args, **kwargs)

        def count_scatter(self, x, y, *args, **kwargs):
            marker = kwargs.get("marker")
            if marker in {"$*$", "$**$", "$***$"}:
                scatter_sizes[marker] = len(x)
            return original_scatter(self, x, y, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "text", fail_star_text)
        monkeypatch.setattr(matplotlib.axes.Axes, "scatter", count_scatter)

        svc._plot_heatmap(corr, pmat, trait_names, out)

        assert scatter_sizes == {"$***$": 2, "$**$": 2, "$*$": 2}
        assert os.path.exists(out)

    def test_draw_significance_stars_no_hits_skips_mask_allocation(
        self, monkeypatch, tmp_path
    ):
        svc = _build_service(str(tmp_path / "corr.png"))
        p_matrix = np.ones((8, 8), dtype=float)

        class Axis:
            def scatter(self, *_args, **_kwargs):
                raise AssertionError("no stars should be drawn")

        def fail_eye(*_args, **_kwargs):
            raise AssertionError("no-star path should not allocate offdiagonal mask")

        def fail_flatnonzero(*_args, **_kwargs):
            raise AssertionError("no-star path should skip coordinate extraction")

        def fail_min(*_args, **_kwargs):
            raise AssertionError("no-star guard should use ndarray.min")

        monkeypatch.setattr(trait_correlation_module.np, "eye", fail_eye)
        monkeypatch.setattr(
            trait_correlation_module.np,
            "min",
            fail_min,
            raising=False,
        )
        monkeypatch.setattr(
            trait_correlation_module.np,
            "flatnonzero",
            fail_flatnonzero,
        )

        svc._draw_significance_stars(Axis(), p_matrix)

    def test_draw_significance_stars_sparse_path_groups_hit_coordinates(
        self, monkeypatch, tmp_path
    ):
        svc = _build_service(str(tmp_path / "corr.png"))
        p_matrix = np.ones((5, 5), dtype=float)
        p_matrix[0, 1] = p_matrix[1, 0] = 0.0005
        p_matrix[0, 2] = p_matrix[2, 0] = 0.005
        p_matrix[3, 4] = p_matrix[4, 3] = 0.02
        p_matrix[2, 2] = 0.0001
        calls = {}

        class Axis:
            def scatter(self, cols, rows, **kwargs):
                calls[kwargs["marker"]] = (cols.tolist(), rows.tolist())

        def fail_eye(*_args, **_kwargs):
            raise AssertionError("sparse star path should not build an eye mask")

        def fail_nonzero(*_args, **_kwargs):
            raise AssertionError("sparse star path should use flat indices")

        monkeypatch.setattr(trait_correlation_module.np, "eye", fail_eye)
        monkeypatch.setattr(trait_correlation_module.np, "nonzero", fail_nonzero)

        svc._draw_significance_stars(Axis(), p_matrix)

        assert calls == {
            "$***$": ([1, 0], [0, 1]),
            "$**$": ([2, 0], [0, 2]),
            "$*$": ([4, 3], [3, 4]),
        }

    def test_json_output(self, tmp_path, capsys):
        out = str(tmp_path / "corr.png")
        args = _make_args(out, json=True)
        svc = TraitCorrelation(args)
        svc.run()

        captured = capsys.readouterr()
        payload = json.loads(captured.out)
        assert "n_taxa" in payload
        assert "n_traits" in payload
        assert "trait_names" in payload
        assert "correlation_matrix" in payload
        assert "p_value_matrix" in payload
        assert "significant_pairs" in payload
        assert "output_file" in payload
        assert payload["n_traits"] == 3
        assert len(payload["correlation_matrix"]) == 3
        assert len(payload["correlation_matrix"][0]) == 3

    def test_print_json_preserves_payload_shape(self, tmp_path, mocker):
        out = str(tmp_path / "corr.png")
        svc = _build_service(out)
        printed_json = mocker.patch(
            "phykit.services.tree.trait_correlation.print_json"
        )
        corr_matrix = np.array(
            [
                [1.0, 0.1234567, -0.5],
                [0.1234567, 1.0, 0.9876543],
                [-0.5, 0.9876543, 1.0],
            ]
        )
        p_matrix = np.array(
            [
                [1.0, 0.0005, 0.02],
                [0.0005, 1.0, 0.2],
                [0.02, 0.2, 1.0],
            ]
        )

        svc._print_json(
            8,
            3,
            ["a", "long_trait", "c"],
            "BM",
            corr_matrix,
            p_matrix,
        )

        printed_json.assert_called_once_with(
            {
                "n_taxa": 8,
                "n_traits": 3,
                "trait_names": ["a", "long_trait", "c"],
                "vcv_type": "BM",
                "alpha": 0.05,
                "correlation_matrix": [
                    [1.0, 0.123457, -0.5],
                    [0.123457, 1.0, 0.987654],
                    [-0.5, 0.987654, 1.0],
                ],
                "p_value_matrix": [
                    [1.0, 0.0005, 0.02],
                    [0.0005, 1.0, 0.2],
                    [0.02, 0.2, 1.0],
                ],
                "significant_pairs": [
                    {
                        "trait_i": "a",
                        "trait_j": "long_trait",
                        "r": 0.123457,
                        "p": 0.0005,
                    },
                    {
                        "trait_i": "a",
                        "trait_j": "c",
                        "r": -0.5,
                        "p": 0.02,
                    },
                ],
                "output_file": out,
            }
        )

    def test_significance_stars(self, tmp_path):
        out = str(tmp_path / "corr.png")
        svc = _build_service(out)
        # Test the stars function
        assert svc._significance_stars(0.0001) == "***"
        assert svc._significance_stars(0.005) == "**"
        assert svc._significance_stars(0.02) == "*"
        assert svc._significance_stars(0.1) == ""

    def test_print_json_reports_upper_triangle_pairs_only(self, tmp_path, mocker):
        out = str(tmp_path / "corr.png")
        svc = _build_service(out)
        printed_json = mocker.patch(
            "phykit.services.tree.trait_correlation.print_json"
        )
        corr_matrix = np.array(
            [
                [1.0, 0.1, 0.2],
                [0.1, 1.0, 0.3],
                [0.2, 0.3, 1.0],
            ]
        )
        p_matrix = np.array(
            [
                [0.001, 0.2, 0.01],
                [0.001, 0.001, 0.2],
                [0.001, 0.001, 0.001],
            ]
        )

        svc._print_json(8, 3, ["a", "b", "c"], "BM", corr_matrix, p_matrix)

        payload = printed_json.call_args.args[0]
        assert payload["significant_pairs"] == [
            {"trait_i": "a", "trait_j": "c", "r": 0.2, "p": 0.01}
        ]

    def test_build_significant_pairs_preserves_upper_triangle_order(self):
        corr_matrix = np.array(
            [
                [1.0, 0.1111111, 0.2222222, 0.3333333],
                [9.0, 1.0, 0.4444444, 0.5555555],
                [9.0, 9.0, 1.0, 0.6666666],
                [9.0, 9.0, 9.0, 1.0],
            ]
        )
        p_matrix = np.array(
            [
                [1.0, 0.01, 0.2, 0.03],
                [0.001, 1.0, 0.04, 0.5],
                [0.001, 0.001, 1.0, 0.02],
                [0.001, 0.001, 0.001, 1.0],
            ]
        )

        pairs = TraitCorrelation._build_significant_pairs(
            ["a", "b", "c", "d"], corr_matrix, p_matrix, 0.05
        )

        assert pairs == [
            {"trait_i": "a", "trait_j": "b", "r": 0.111111, "p": 0.01},
            {"trait_i": "a", "trait_j": "d", "r": 0.333333, "p": 0.03},
            {"trait_i": "b", "trait_j": "c", "r": 0.444444, "p": 0.04},
            {"trait_i": "c", "trait_j": "d", "r": 0.666667, "p": 0.02},
        ]

    def test_build_significant_pairs_sparse_path_avoids_triangular_mask(
        self, monkeypatch
    ):
        trait_names = [f"trait_{idx}" for idx in range(20)]
        corr_matrix = np.zeros((20, 20), dtype=float)
        p_matrix = np.ones((20, 20), dtype=float)
        corr_matrix[2, 19] = 0.1234567
        p_matrix[2, 19] = 0.01

        def fail_triu(*_args, **_kwargs):
            raise AssertionError("sparse significant pairs should scan rows directly")

        monkeypatch.setattr(trait_correlation_module.np, "triu", fail_triu)

        pairs = TraitCorrelation._build_significant_pairs(
            trait_names, corr_matrix, p_matrix, 0.05
        )

        assert pairs == [
            {
                "trait_i": "trait_2",
                "trait_j": "trait_19",
                "r": 0.123457,
                "p": 0.01,
            }
        ]

    def test_text_output_format(self, tmp_path, capsys):
        out = str(tmp_path / "corr.png")
        args = _make_args(out, json=False)
        svc = TraitCorrelation(args)
        svc.run()

        captured = capsys.readouterr()
        text = captured.out
        assert "Phylogenetic Trait Correlation" in text
        assert "body_mass" in text
        assert "brain_size" in text
        assert "longevity" in text
        assert "Taxa:" in text
        assert "Traits:" in text
        assert "Significant pairs" in text

    def test_print_text_batches_matrix_rows_and_counts_upper_triangle(
        self, tmp_path, mocker
    ):
        out = str(tmp_path / "corr.png")
        svc = _build_service(out)
        printed = mocker.patch("builtins.print")
        corr_matrix = np.array(
            [
                [1.0, 0.12345, -0.5],
                [0.12345, 1.0, 0.98765],
                [-0.5, 0.98765, 1.0],
            ]
        )
        p_matrix = np.array(
            [
                [1.0, 0.0005, 0.02],
                [0.0005, 1.0, 0.2],
                [0.02, 0.2, 1.0],
            ]
        )
        mocker.patch.object(
            svc,
            "_significance_stars",
            side_effect=AssertionError(
                "_print_text should inline star thresholds while formatting rows"
            ),
        )
        mocker.patch.object(
            trait_correlation_module.np,
            "triu_indices",
            side_effect=AssertionError(
                "_print_text should count significant pairs without triangular index arrays"
            ),
        )

        svc._print_text(
            8,
            3,
            ["a", "long_trait", "c"],
            "BM",
            corr_matrix,
            p_matrix,
            "tree.nwk",
        )

        printed.assert_called_once_with(
            "Phylogenetic Trait Correlation\n"
            "Tree: tree.nwk\n"
            "Taxa: 8\n"
            "Traits: 3\n"
            "VCV: BM (standard)\n"
            "Alpha: 0.05\n"
            "\n"
            "                       a    long_trait             c  \n"
            "a                 1.0000     0.1235***      -0.5000*  \n"
            "long_trait     0.1235***        1.0000        0.9877  \n"
            "c               -0.5000*        0.9877        1.0000  \n"
            "\n"
            "Significant pairs (p < 0.05): 2 of 3\n"
            f"Output: {svc.output_path}"
        )

    def test_count_significant_upper_triangle_scans_rows(self, monkeypatch):
        p_matrix = np.array(
            [
                [1.0, 0.01, 0.2, 0.03],
                [0.001, 1.0, 0.04, 0.5],
                [0.001, 0.001, 1.0, 0.02],
                [0.001, 0.001, 0.001, 1.0],
            ]
        )

        def fail_triu_indices(*_args, **_kwargs):
            raise AssertionError("significance count should scan upper-triangle rows")

        monkeypatch.setattr(
            trait_correlation_module.np,
            "triu_indices",
            fail_triu_indices,
        )

        assert TraitCorrelation._count_significant_upper_triangle(p_matrix, 0.05) == 4

    def _compute_correlation(self, svc):
        """Helper: re-compute the correlation matrix from a built service."""
        from phykit.services.tree.vcv_utils import build_vcv_matrix
        from scipy.stats import t as t_dist

        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(
            svc.trait_data_path, tree_tips, min_columns=3
        )
        ordered_names = sorted(traits.keys())
        vcv = build_vcv_matrix(tree, ordered_names)

        n = len(ordered_names)
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])

        C_inv = np.linalg.inv(vcv)
        ones = np.ones(n)
        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])
        Y_centered = Y - a_hat

        phylo_cov = (Y_centered.T @ C_inv @ Y_centered) / (n - 1)
        std_devs = np.sqrt(np.diag(phylo_cov))
        corr_matrix = phylo_cov / np.outer(std_devs, std_devs)
        np.fill_diagonal(corr_matrix, 1.0)

        p_matrix = np.ones((p, p))
        for i in range(p):
            for j in range(i + 1, p):
                r = corr_matrix[i, j]
                if abs(r) >= 1.0 - 1e-10:
                    pval = 0.0
                else:
                    t_stat = r * np.sqrt((n - 2) / (1 - r ** 2))
                    pval = 2.0 * t_dist.sf(abs(t_stat), df=n - 2)
                p_matrix[i, j] = pval
                p_matrix[j, i] = pval

        return corr_matrix, p_matrix, trait_names, n
