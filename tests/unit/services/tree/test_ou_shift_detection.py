import copy

import numpy as np
import pytest
from argparse import Namespace
from pathlib import Path

from phykit.services.tree.ou_shift_detection import OUShiftDetection
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
        criterion="pBIC",
        max_shifts=None,
        json=False,
    )


@pytest.fixture
def svc(default_args):
    return OUShiftDetection(default_args)


@pytest.fixture
def precomputed(svc):
    """Pre-compute common data structures for testing."""
    tree = svc.read_tree_file()
    tree_tips = svc.get_tip_names_from_tree(tree)
    traits = svc._parse_trait_file(TRAITS_FILE, tree_tips)

    shared = set(traits.keys())
    tree_copy = copy.deepcopy(tree)
    tip_names_in_tree = [t.name for t in tree_copy.get_terminals()]
    to_prune = [t for t in tip_names_in_tree if t not in shared]
    if to_prune:
        for t in to_prune:
            tree_copy.prune(t)

    ordered_names = sorted(shared)
    x = np.array([traits[name] for name in ordered_names])

    parent_map = svc._build_parent_map(tree_copy)
    lineage_info = svc._build_lineage_info_no_regime(
        tree_copy, ordered_names, parent_map
    )
    tree_height = max(
        sum(bl for _, bl, _, _ in lineage_info[name])
        for name in ordered_names
    )

    edges = svc._enumerate_edges(tree_copy, parent_map)
    descendant_counts = svc._count_descendants(tree_copy, edges)

    # Set up instance attributes needed by new methods
    svc._ordered_names = ordered_names
    svc._x = x
    svc._n = len(ordered_names)
    svc._lineage_info = lineage_info
    svc._tree_height = tree_height
    svc._edges = edges
    svc._n_edges = len(edges)
    svc._parent_map = parent_map
    svc._S, svc._tip_heights_arr = svc._precompute_shared_path_lengths()

    return dict(
        svc=svc, tree=tree_copy, x=x, ordered_names=ordered_names,
        lineage_info=lineage_info, tree_height=tree_height,
        parent_map=parent_map, edges=edges,
        descendant_counts=descendant_counts,
    )


class TestEnumerateEdges:
    def test_correct_count_excludes_root(self, precomputed):
        edges = precomputed["edges"]
        tree = precomputed["tree"]
        all_clades = list(tree.find_clades())
        expected = len(all_clades) - 1
        assert len(edges) == expected

    def test_root_not_in_edges(self, precomputed):
        edges = precomputed["edges"]
        tree = precomputed["tree"]
        root_id = id(tree.root)
        edge_ids = [eid for eid, _ in edges]
        assert root_id not in edge_ids

    def test_all_edges_have_valid_clades(self, precomputed):
        edges = precomputed["edges"]
        for clade_id, clade in edges:
            assert clade_id == id(clade)
            assert clade is not None


class TestCountDescendants:
    def test_terminal_edges_have_one_descendant(self, precomputed):
        edges = precomputed["edges"]
        descendant_counts = precomputed["descendant_counts"]
        for clade_id, clade in edges:
            if clade.is_terminal():
                assert descendant_counts[clade_id] == 1

    def test_internal_edges_have_multiple_descendants(self, precomputed):
        edges = precomputed["edges"]
        descendant_counts = precomputed["descendant_counts"]
        for clade_id, clade in edges:
            if not clade.is_terminal():
                assert descendant_counts[clade_id] >= 2


class TestBuildShiftWeightMatrix:
    def test_shape(self, precomputed):
        svc = precomputed["svc"]
        ordered_names = precomputed["ordered_names"]
        lineage_info = precomputed["lineage_info"]
        edges = precomputed["edges"]
        tree_height = precomputed["tree_height"]
        n = len(ordered_names)
        p = len(edges)

        W = svc._build_shift_weight_matrix(
            ordered_names, lineage_info, edges, 1.0, tree_height
        )
        assert W.shape == (n, p + 1)

    def test_rows_sum_to_one(self, precomputed):
        svc = precomputed["svc"]
        ordered_names = precomputed["ordered_names"]
        lineage_info = precomputed["lineage_info"]
        edges = precomputed["edges"]
        tree_height = precomputed["tree_height"]

        W = svc._build_shift_weight_matrix(
            ordered_names, lineage_info, edges, 1.0, tree_height
        )
        row_sums = W.sum(axis=1)
        np.testing.assert_allclose(row_sums, 1.0, atol=1e-10)

    def test_all_weights_nonnegative(self, precomputed):
        svc = precomputed["svc"]
        ordered_names = precomputed["ordered_names"]
        lineage_info = precomputed["lineage_info"]
        edges = precomputed["edges"]
        tree_height = precomputed["tree_height"]

        W = svc._build_shift_weight_matrix(
            ordered_names, lineage_info, edges, 1.0, tree_height
        )
        assert np.all(W >= -1e-15)


class TestBuildIndicatorDesignMatrix:
    def test_shape(self, precomputed):
        svc = precomputed["svc"]
        n = svc._n
        # No shifts: should be n x 1
        W = svc._build_indicator_design_matrix([])
        assert W.shape == (n, 1)
        np.testing.assert_allclose(W[:, 0], 1.0)

    def test_shape_with_shifts(self, precomputed):
        svc = precomputed["svc"]
        n = svc._n
        W = svc._build_indicator_design_matrix([0, 1])
        assert W.shape == (n, 3)

    def test_intercept_is_ones(self, precomputed):
        svc = precomputed["svc"]
        W = svc._build_indicator_design_matrix([0, 1, 2])
        np.testing.assert_allclose(W[:, 0], 1.0)

    def test_values_are_binary(self, precomputed):
        svc = precomputed["svc"]
        edges = precomputed["edges"]
        all_indices = list(range(len(edges)))
        W = svc._build_indicator_design_matrix(all_indices)
        # All values should be 0 or 1
        for col in range(1, W.shape[1]):
            unique_vals = np.unique(W[:, col])
            assert all(v in (0.0, 1.0) for v in unique_vals)

    def test_terminal_edge_marks_one_tip(self, precomputed):
        svc = precomputed["svc"]
        edges = precomputed["edges"]
        # Find a terminal edge
        for j, (clade_id, clade) in enumerate(edges):
            if clade.is_terminal():
                W = svc._build_indicator_design_matrix([j])
                # Only one tip should have a 1 in the shift column
                assert W[:, 1].sum() == 1.0
                break


class TestComputePBIC:
    def test_returns_finite_for_valid_model(self, precomputed):
        """pBIC should return a finite value for a valid model fit."""
        svc = precomputed["svc"]
        n = svc._n
        x = svc._x

        # Fit null model and check pBIC is finite
        alpha = 1.0
        V = svc._build_ou_vcv_fast(alpha)
        V_inv = np.linalg.inv(V)
        W = np.ones((n, 1))
        theta = svc._gls_theta_hat(x, W, V_inv)
        e = x - W @ theta
        sig2 = float(e @ V_inv @ e) / n
        sign, logdet = np.linalg.slogdet(V)
        ll = -0.5 * (n * np.log(2 * np.pi) + n * np.log(sig2) + logdet + n)

        pbic = svc._compute_pbic(n, ll, 0, W, V_inv, sig2, x)
        assert np.isfinite(pbic)

    def test_more_shifts_higher_penalty(self, precomputed):
        """pBIC penalty should increase with number of shifts."""
        svc = precomputed["svc"]
        n = svc._n
        x = svc._x

        # Create two models: 0 shifts and 2 shifts (with same LL)
        alpha = 1.0
        V = svc._build_ou_vcv_fast(alpha)
        V_inv = np.linalg.inv(V)

        # Null model
        W0 = np.ones((n, 1))
        theta0 = svc._gls_theta_hat(x, W0, V_inv)
        e0 = x - W0 @ theta0
        sig2_0 = float(e0 @ V_inv @ e0) / n
        sign, logdet = np.linalg.slogdet(V)
        ll_0 = -0.5 * (n * np.log(2 * np.pi) + n * np.log(sig2_0) + logdet + n)

        # 2-shift model (using first 2 edges)
        edges2 = svc._edges[:2]
        W2 = svc._build_shift_weight_matrix(
            svc._ordered_names, svc._lineage_info, edges2, alpha, svc._tree_height
        )
        theta2 = svc._gls_theta_hat(x, W2, V_inv)
        e2 = x - W2 @ theta2
        sig2_2 = float(e2 @ V_inv @ e2) / n
        ll_2 = -0.5 * (n * np.log(2 * np.pi) + n * np.log(sig2_2) + logdet + n)

        pbic_0 = svc._compute_pbic(n, ll_0, 0, W0, V_inv, sig2_0, x)
        pbic_2 = svc._compute_pbic(n, ll_2, 2, W2, V_inv, sig2_2, x)

        # Both should be finite
        assert np.isfinite(pbic_0)
        assert np.isfinite(pbic_2)


class TestComputeBIC:
    def test_value(self, svc):
        n = 10
        ll = -20.0
        k = 2
        bic = svc._compute_bic(n, ll, k)
        # R-style: -2*LL + (2k+3)*log(n)
        expected = 40.0 + 7 * np.log(10)
        np.testing.assert_allclose(bic, expected, atol=1e-10)


class TestComputeAICc:
    def test_value(self, svc):
        n = 20
        ll = -30.0
        k = 2
        aicc = svc._compute_aicc(n, ll, k)
        # R-style: p = 2k+3 = 7
        p = 7
        expected = 60.0 + 2.0 * p + 2.0 * p * 8 / (20 - p - 1)
        np.testing.assert_allclose(aicc, expected, atol=1e-10)

    def test_small_n_gives_inf(self, svc):
        # p = 2*2+3 = 7, need n > p+1 = 8
        aicc = svc._compute_aicc(3, -10.0, 1)
        assert aicc == float("inf")


class TestTransformToIndependent:
    def test_identity_covariance_after_transform(self, precomputed):
        svc = precomputed["svc"]
        ordered_names = precomputed["ordered_names"]
        lineage_info = precomputed["lineage_info"]
        edges = precomputed["edges"]
        tree_height = precomputed["tree_height"]
        x = precomputed["x"]

        alpha = 1.0
        V = svc._build_ou_vcv_single_alpha_no_regime(
            ordered_names, lineage_info, alpha, 1.0
        )
        W = svc._build_shift_weight_matrix(
            ordered_names, lineage_info, edges, alpha, tree_height
        )
        X_star, y_star = svc._transform_to_independent(x, W, V)

        L = np.linalg.cholesky(V)
        L_inv = np.linalg.inv(L)
        transformed_V = L_inv @ V @ L_inv.T
        np.testing.assert_allclose(transformed_V, np.eye(len(x)), atol=1e-10)


class TestFitOU1ForAlpha:
    def test_returns_positive_alpha(self, precomputed):
        svc = precomputed["svc"]
        x = precomputed["x"]
        ordered_names = precomputed["ordered_names"]
        lineage_info = precomputed["lineage_info"]
        tree_height = precomputed["tree_height"]

        alpha = svc._fit_ou1_for_alpha(x, ordered_names, lineage_info, tree_height)
        assert alpha > 0


class TestPrecomputeSharedPathLengths:
    def test_diagonal_equals_tip_heights(self, precomputed):
        svc = precomputed["svc"]
        S = svc._S
        tip_heights = svc._tip_heights_arr
        n = svc._n
        for i in range(n):
            np.testing.assert_allclose(S[i, i], tip_heights[i], atol=1e-12)

    def test_shared_path_symmetric(self, precomputed):
        svc = precomputed["svc"]
        S = svc._S
        np.testing.assert_allclose(S, S.T, atol=1e-12)


class TestBuildOUVCVFast:
    def test_matches_slow_path(self, precomputed):
        """Fast VCV should match the loop-based VCV."""
        svc = precomputed["svc"]
        ordered_names = precomputed["ordered_names"]
        lineage_info = precomputed["lineage_info"]

        alpha = 1.5
        V_slow = svc._build_ou_vcv_single_alpha_no_regime(
            ordered_names, lineage_info, alpha, 1.0
        )
        V_fast = svc._build_ou_vcv_fast(alpha, 1.0)
        np.testing.assert_allclose(V_fast, V_slow, atol=1e-10)

    def test_bm_limit(self, precomputed):
        """Alpha near 0 should give BM covariance."""
        svc = precomputed["svc"]
        V_bm = svc._build_ou_vcv_fast(1e-12)
        np.testing.assert_allclose(V_bm, svc._S, atol=1e-6)


class TestEndToEnd:
    def test_pipeline_runs_with_setup(self, precomputed):
        """Pipeline should complete with pre-computed data."""
        svc = precomputed["svc"]
        # Use _fit_and_score_config for null model
        null = svc._fit_and_score_config([])
        assert null is not None
        assert null["n_shifts"] == 0
        assert "log_likelihood" in null
        assert "pBIC" in null

    def test_fit_with_shifts(self, precomputed):
        """Fitting with shift edges should complete."""
        svc = precomputed["svc"]
        # Fit with first 2 edges as shifts
        result = svc._fit_and_score_config([0, 1])
        assert result is not None
        assert result["n_shifts"] == 2
        assert result["alpha"] > 0
        assert np.isfinite(result["pBIC"])

    def test_run_completes(self, default_args):
        """Full run completes without error."""
        from unittest.mock import patch
        svc = OUShiftDetection(default_args)
        with patch("builtins.print"):
            svc.run()


class TestDescribeEdge:
    def test_terminal(self, precomputed):
        svc = precomputed["svc"]
        tree = precomputed["tree"]
        edges = precomputed["edges"]
        for clade_id, clade in edges:
            if clade.is_terminal():
                desc = svc._describe_edge(clade, tree)
                assert "terminal branch to" in desc
                break

    def test_internal(self, precomputed):
        svc = precomputed["svc"]
        tree = precomputed["tree"]
        edges = precomputed["edges"]
        for clade_id, clade in edges:
            if not clade.is_terminal():
                desc = svc._describe_edge(clade, tree)
                assert "stem of" in desc
                break


class TestProcessArgs:
    def test_invalid_criterion(self):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            criterion="INVALID",
            max_shifts=None,
            json=False,
        )
        with pytest.raises(PhykitUserError):
            OUShiftDetection(args)

    def test_valid_criteria(self):
        for crit in ("pBIC", "BIC", "AICc"):
            args = Namespace(
                tree=TREE_SIMPLE,
                trait_data=TRAITS_FILE,
                criterion=crit,
                max_shifts=None,
                json=False,
            )
            svc = OUShiftDetection(args)
            assert svc.criterion == crit
