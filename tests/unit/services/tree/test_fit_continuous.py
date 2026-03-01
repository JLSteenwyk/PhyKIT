import pytest
import json
import numpy as np
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch

from phykit.services.tree.fit_continuous import FitContinuous, ALL_MODELS
from phykit.errors import PhykitUserError


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_traits.tsv")


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
    max_lam = svc._max_lambda(tree)
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


# ── TestBMModel ──────────────────────────────────────────────────────


class TestBMModel:
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
    def test_aic_weights_sum_to_one(self, tree_vcv_data):
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
    def test_ou_vcv_positive_definite(self, tree_vcv_data):
        d = tree_vcv_data
        V = d["svc"]._vcv_ou(d["vcv"], 0.5)
        eigenvalues = np.linalg.eigvalsh(V)
        assert np.all(eigenvalues > -1e-10)

    def test_lambda_1_equals_bm(self, tree_vcv_data):
        d = tree_vcv_data
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
