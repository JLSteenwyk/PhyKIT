import json

import numpy as np
import pytest
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch

from phykit.services.tree.ouwie import OUwie, ALL_MODELS
from phykit.errors import PhykitUserError


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_traits.tsv")
REGIMES_FILE = str(SAMPLE_FILES / "tree_simple_regimes.tsv")


@pytest.fixture
def default_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=TRAITS_FILE,
        regime_data=REGIMES_FILE,
        models=None,
        json=False,
    )


@pytest.fixture
def svc(default_args):
    return OUwie(default_args)


@pytest.fixture
def precomputed(svc):
    """Pre-compute all common data structures."""
    import copy
    tree = svc.read_tree_file()
    tree_tips = svc.get_tip_names_from_tree(tree)
    traits = svc._parse_trait_file(TRAITS_FILE, tree_tips)
    regime_assignments = svc._parse_regime_file(REGIMES_FILE, tree_tips)

    shared = set(traits.keys()) & set(regime_assignments.keys())
    traits = {k: traits[k] for k in shared}
    regime_assignments = {k: regime_assignments[k] for k in shared}

    tree_copy = copy.deepcopy(tree)
    tip_names_in_tree = [t.name for t in tree_copy.get_terminals()]
    tips_to_prune = [t for t in tip_names_in_tree if t not in shared]
    if tips_to_prune:
        tree_copy = svc.prune_tree_using_taxa_list(tree_copy, tips_to_prune)

    ordered_names = sorted(traits.keys())
    x = np.array([traits[name] for name in ordered_names])

    regimes = sorted(set(regime_assignments.values()))
    parent_map = svc._build_parent_map(tree_copy)
    branch_regimes = svc._assign_branch_regimes(
        tree_copy, regime_assignments, parent_map
    )
    root_regime = svc._get_root_regime(tree_copy, regime_assignments, parent_map)

    lineage_info = svc._build_lineage_info(
        tree_copy, ordered_names, parent_map, branch_regimes
    )

    per_regime_vcv = svc._build_per_regime_vcv(
        tree_copy, ordered_names, regimes, branch_regimes, parent_map
    )
    vcv_total = sum(per_regime_vcv.values())
    tree_height = float(np.max(np.diag(vcv_total)))

    return dict(
        svc=svc, tree=tree_copy, x=x, ordered_names=ordered_names,
        regimes=regimes, lineage_info=lineage_info,
        per_regime_vcv=per_regime_vcv, vcv_total=vcv_total,
        root_regime=root_regime, tree_height=tree_height,
        n_regimes=len(regimes),
    )


# ── TestProcessArgs ──────────────────────────────────────────────────

class TestProcessArgs:
    def test_defaults(self):
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=TRAITS_FILE,
            regime_data=REGIMES_FILE, models=None, json=False,
        )
        svc = OUwie(args)
        assert svc.selected_models == ALL_MODELS

    def test_subset_models(self):
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=TRAITS_FILE,
            regime_data=REGIMES_FILE, models="BM1,OUM,OUMVA", json=False,
        )
        svc = OUwie(args)
        assert svc.selected_models == ["BM1", "OUM", "OUMVA"]

    def test_invalid_model(self):
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=TRAITS_FILE,
            regime_data=REGIMES_FILE, models="BM1,BOGUS", json=False,
        )
        with pytest.raises(PhykitUserError):
            OUwie(args)


# ── TestWeightMatrix ─────────────────────────────────────────────────

class TestWeightMatrix:
    def test_rows_sum_to_one(self, precomputed):
        d = precomputed
        alpha = 0.5
        W = d["svc"]._build_weight_matrix_single_alpha(
            d["ordered_names"], d["lineage_info"], d["regimes"],
            alpha, d["root_regime"],
        )
        row_sums = W.sum(axis=1)
        np.testing.assert_allclose(row_sums, 1.0, atol=1e-10)

    def test_single_regime_reduces_to_ones(self, precomputed):
        d = precomputed
        alpha = 1.0
        # If all taxa were same regime, W would be n x 1 of ones
        # With 2 regimes, rows should still sum to 1
        W = d["svc"]._build_weight_matrix_single_alpha(
            d["ordered_names"], d["lineage_info"], d["regimes"],
            alpha, d["root_regime"],
        )
        assert W.shape == (len(d["ordered_names"]), len(d["regimes"]))
        assert np.all(W >= -1e-10)  # Non-negative weights

    def test_multi_alpha_rows_sum_to_one(self, precomputed):
        d = precomputed
        alphas_dict = {r: 0.5 for r in d["regimes"]}
        W = d["svc"]._build_weight_matrix_multi_alpha(
            d["ordered_names"], d["lineage_info"], d["regimes"],
            alphas_dict, d["root_regime"],
        )
        row_sums = W.sum(axis=1)
        np.testing.assert_allclose(row_sums, 1.0, atol=1e-10)


# ── TestVCV ──────────────────────────────────────────────────────────

class TestVCV:
    def test_positive_definite(self, precomputed):
        d = precomputed
        alpha = 0.5
        V = d["svc"]._build_ou_vcv_single_alpha(
            d["ordered_names"], d["lineage_info"], alpha, 1.0,
        )
        eigenvalues = np.linalg.eigvalsh(V)
        assert np.all(eigenvalues > 0)

    def test_H_matrices_sum_correctly(self, precomputed):
        d = precomputed
        alpha = 0.5
        H = d["svc"]._build_ou_H_matrices(
            d["ordered_names"], d["lineage_info"], d["regimes"], alpha,
        )
        # Sum of H_r * sigma2=1 should equal single-alpha VCV with sigma2=1
        V_from_H = sum(H[r] for r in d["regimes"])
        V_direct = d["svc"]._build_ou_vcv_single_alpha(
            d["ordered_names"], d["lineage_info"], alpha, 1.0,
        )
        np.testing.assert_allclose(V_from_H, V_direct, atol=1e-10)

    def test_multi_alpha_positive_definite(self, precomputed):
        d = precomputed
        alphas_dict = {r: 0.5 for r in d["regimes"]}
        sigmas_dict = {r: 1.0 for r in d["regimes"]}
        V = d["svc"]._build_ou_vcv_multi_alpha(
            d["ordered_names"], d["lineage_info"], d["regimes"],
            alphas_dict, sigmas_dict,
        )
        eigenvalues = np.linalg.eigvalsh(V)
        assert np.all(eigenvalues > 0)


# ── TestBM1 ──────────────────────────────────────────────────────────

class TestBM1:
    def test_k_params_equals_2(self, precomputed):
        d = precomputed
        res = d["svc"]._fit_bm1(d["x"], d["vcv_total"])
        assert res["k_params"] == 2

    def test_ll_is_finite(self, precomputed):
        d = precomputed
        res = d["svc"]._fit_bm1(d["x"], d["vcv_total"])
        assert np.isfinite(res["log_likelihood"])


# ── TestBMS ──────────────────────────────────────────────────────────

class TestBMS:
    def test_ll_ge_bm1(self, precomputed):
        d = precomputed
        bm1 = d["svc"]._fit_bm1(d["x"], d["vcv_total"])
        bms = d["svc"]._fit_bms(d["x"], d["per_regime_vcv"], d["regimes"])
        # BMS has more parameters, should fit at least as well
        assert bms["log_likelihood"] >= bm1["log_likelihood"] - 1e-6

    def test_per_regime_sigma2_positive(self, precomputed):
        d = precomputed
        bms = d["svc"]._fit_bms(d["x"], d["per_regime_vcv"], d["regimes"])
        for r in d["regimes"]:
            assert bms["params"]["sigma2"][r] > 0


# ── TestOU1 ──────────────────────────────────────────────────────────

class TestOU1:
    def test_alpha_positive(self, precomputed):
        d = precomputed
        res = d["svc"]._fit_ou1(
            d["x"], d["ordered_names"], d["lineage_info"], d["tree_height"],
        )
        assert res["params"]["alpha"] > 0

    def test_ll_is_finite(self, precomputed):
        d = precomputed
        res = d["svc"]._fit_ou1(
            d["x"], d["ordered_names"], d["lineage_info"], d["tree_height"],
        )
        assert np.isfinite(res["log_likelihood"])


# ── TestOUM ──────────────────────────────────────────────────────────

class TestOUM:
    def test_ll_ge_ou1(self, precomputed):
        d = precomputed
        ou1 = d["svc"]._fit_ou1(
            d["x"], d["ordered_names"], d["lineage_info"], d["tree_height"],
        )
        oum = d["svc"]._fit_oum(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        assert oum["log_likelihood"] >= ou1["log_likelihood"] - 1e-6

    def test_per_regime_theta(self, precomputed):
        d = precomputed
        oum = d["svc"]._fit_oum(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        assert isinstance(oum["params"]["theta"], dict)
        assert len(oum["params"]["theta"]) == len(d["regimes"])


# ── TestOUMV ─────────────────────────────────────────────────────────

class TestOUMV:
    def test_ll_ge_oum(self, precomputed):
        d = precomputed
        oum = d["svc"]._fit_oum(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        oumv = d["svc"]._fit_oumv(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        assert oumv["log_likelihood"] >= oum["log_likelihood"] - 1e-6

    def test_per_regime_sigma2(self, precomputed):
        d = precomputed
        oumv = d["svc"]._fit_oumv(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        assert isinstance(oumv["params"]["sigma2"], dict)
        for r in d["regimes"]:
            assert oumv["params"]["sigma2"][r] > 0


# ── TestOUMA ─────────────────────────────────────────────────────────

class TestOUMA:
    def test_ll_ge_oum(self, precomputed):
        d = precomputed
        oum = d["svc"]._fit_oum(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        ouma = d["svc"]._fit_ouma(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        assert ouma["log_likelihood"] >= oum["log_likelihood"] - 1e-6

    def test_per_regime_alpha(self, precomputed):
        d = precomputed
        ouma = d["svc"]._fit_ouma(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        assert isinstance(ouma["params"]["alpha"], dict)
        for r in d["regimes"]:
            assert ouma["params"]["alpha"][r] > 0


# ── TestOUMVA ────────────────────────────────────────────────────────

class TestOUMVA:
    def test_ll_ge_oumv_and_ouma(self, precomputed):
        d = precomputed
        oumv = d["svc"]._fit_oumv(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        ouma = d["svc"]._fit_ouma(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        oumva = d["svc"]._fit_oumva(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        assert oumva["log_likelihood"] >= oumv["log_likelihood"] - 1e-6
        assert oumva["log_likelihood"] >= ouma["log_likelihood"] - 1e-6

    def test_k_params_equals_3R(self, precomputed):
        d = precomputed
        R = len(d["regimes"])
        oumva = d["svc"]._fit_oumva(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        assert oumva["k_params"] == 3 * R


# ── TestModelNesting ─────────────────────────────────────────────────

class TestModelNesting:
    """Verify that nested models have LL(simpler) <= LL(complex)."""

    def test_bm1_le_bms(self, precomputed):
        d = precomputed
        bm1 = d["svc"]._fit_bm1(d["x"], d["vcv_total"])
        bms = d["svc"]._fit_bms(d["x"], d["per_regime_vcv"], d["regimes"])
        assert bms["log_likelihood"] >= bm1["log_likelihood"] - 1e-6

    def test_ou1_le_oum(self, precomputed):
        d = precomputed
        ou1 = d["svc"]._fit_ou1(
            d["x"], d["ordered_names"], d["lineage_info"], d["tree_height"],
        )
        oum = d["svc"]._fit_oum(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        assert oum["log_likelihood"] >= ou1["log_likelihood"] - 1e-6

    def test_oum_le_oumv(self, precomputed):
        d = precomputed
        oum = d["svc"]._fit_oum(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        oumv = d["svc"]._fit_oumv(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        assert oumv["log_likelihood"] >= oum["log_likelihood"] - 1e-6

    def test_oum_le_ouma(self, precomputed):
        d = precomputed
        oum = d["svc"]._fit_oum(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        ouma = d["svc"]._fit_ouma(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        assert ouma["log_likelihood"] >= oum["log_likelihood"] - 1e-6

    def test_oumv_le_oumva(self, precomputed):
        d = precomputed
        oumv = d["svc"]._fit_oumv(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        oumva = d["svc"]._fit_oumva(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        assert oumva["log_likelihood"] >= oumv["log_likelihood"] - 1e-6

    def test_ouma_le_oumva(self, precomputed):
        d = precomputed
        ouma = d["svc"]._fit_ouma(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        oumva = d["svc"]._fit_oumva(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        assert oumva["log_likelihood"] >= ouma["log_likelihood"] - 1e-6


# ── TestModelComparison ──────────────────────────────────────────────

class TestModelComparison:
    def test_aicc_weights_sum_to_one(self, precomputed):
        d = precomputed
        results = []
        results.append(d["svc"]._fit_bm1(d["x"], d["vcv_total"]))
        results.append(d["svc"]._fit_ou1(
            d["x"], d["ordered_names"], d["lineage_info"], d["tree_height"],
        ))
        results.append(d["svc"]._fit_oum(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        ))
        n = len(d["x"])
        results = d["svc"]._compute_model_comparison(results, n)
        weights = [r["aicc_weight"] for r in results]
        np.testing.assert_allclose(sum(weights), 1.0, atol=1e-10)

    def test_sorted_by_aicc(self, precomputed):
        d = precomputed
        results = []
        results.append(d["svc"]._fit_bm1(d["x"], d["vcv_total"]))
        results.append(d["svc"]._fit_ou1(
            d["x"], d["ordered_names"], d["lineage_info"], d["tree_height"],
        ))
        n = len(d["x"])
        results = d["svc"]._compute_model_comparison(results, n)
        aiccs = [r["aicc"] for r in results]
        assert aiccs == sorted(aiccs)

    def test_all_lls_finite(self, precomputed):
        d = precomputed
        results = []
        for model in ALL_MODELS:
            res = d["svc"]._fit_model(
                model, d["x"], d["vcv_total"], d["per_regime_vcv"],
                d["ordered_names"], d["regimes"], d["lineage_info"],
                d["root_regime"], d["tree_height"], d["n_regimes"],
            )
            results.append(res)
        for r in results:
            assert np.isfinite(r["log_likelihood"]), f"{r['model']} LL not finite"


# ── TestRun ──────────────────────────────────────────────────────────

class TestRun:
    @patch("builtins.print")
    def test_text_output(self, mocked_print, svc):
        svc.run()
        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "OUwie Model Comparison" in all_output
        assert "BM1" in all_output
        assert "OUM" in all_output

    @patch("builtins.print")
    def test_json_output(self, mocked_print):
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=TRAITS_FILE,
            regime_data=REGIMES_FILE, models=None, json=True,
        )
        svc = OUwie(args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["n_tips"] == 8
        assert "models" in payload
        assert "BM1" in payload["models"]
        assert "best_model_aicc" in payload

    @patch("builtins.print")
    def test_subset_models(self, mocked_print):
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=TRAITS_FILE,
            regime_data=REGIMES_FILE, models="BM1,OUM", json=True,
        )
        svc = OUwie(args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        assert set(payload["models"].keys()) == {"BM1", "OUM"}
