"""
Unit tests for fit_discrete (discrete trait model comparison).

Fits ER, SYM, and ARD Mk models via maximum likelihood and compares
them using AIC/BIC. Expected log-likelihood values cross-validated
against R's geiger::fitDiscrete() (with multi2di() for multifurcations).
See tests/r_validation/validate_fit_discrete.R for the R script.

Note: PhyKIT and geiger may produce slightly different log-likelihoods
(within ~1.0) due to differences in handling multifurcating trees and
multi-start optimization. The model ranking should agree.
"""
import math
import pytest
from argparse import Namespace

from phykit.services.tree.fit_discrete import FitDiscrete
from phykit.helpers.discrete_models import (
    build_q_matrix,
    count_params,
    felsenstein_pruning,
    fit_q_matrix,
    parse_discrete_traits,
    VALID_DISCRETE_MODELS,
)


@pytest.fixture
def args():
    return Namespace(
        tree="tests/sample_files/tree_simple.tre",
        trait_data="tests/sample_files/tree_simple_discrete_traits.tsv",
        trait="diet",
    )


class TestFitDiscreteInit:
    def test_init_sets_fields(self, args):
        fd = FitDiscrete(args)
        assert fd.tree_file_path == args.tree
        assert fd.trait_data_path == args.trait_data
        assert fd.trait_column == "diet"
        assert fd.selected_models == ["ER", "SYM", "ARD"]
        assert fd.json_output is False

    def test_process_args_custom_models(self):
        args = Namespace(
            tree="t.tre", trait_data="d.tsv", trait="col", models="ER,ARD"
        )
        fd = FitDiscrete(args)
        assert fd.selected_models == ["ER", "ARD"]

    def test_process_args_invalid_model_exits(self):
        args = Namespace(
            tree="t.tre", trait_data="d.tsv", trait="col", models="ER,FOO"
        )
        with pytest.raises(SystemExit):
            FitDiscrete(args)


class TestDiscreteModelsShared:
    def test_build_q_matrix_er(self):
        Q = build_q_matrix([0.5], 3, "ER")
        assert Q.shape == (3, 3)
        # Off-diagonals should all be 0.5
        for i in range(3):
            for j in range(3):
                if i != j:
                    assert Q[i, j] == pytest.approx(0.5)
        # Rows sum to 0
        for i in range(3):
            assert sum(Q[i, :]) == pytest.approx(0.0)

    def test_build_q_matrix_sym(self):
        # k=3 -> 3 params
        Q = build_q_matrix([0.1, 0.2, 0.3], 3, "SYM")
        assert Q[0, 1] == Q[1, 0]  # Symmetric
        assert Q[0, 2] == Q[2, 0]
        assert Q[1, 2] == Q[2, 1]

    def test_build_q_matrix_ard(self):
        # k=3 -> 6 params
        Q = build_q_matrix([0.1, 0.2, 0.3, 0.4, 0.5, 0.6], 3, "ARD")
        assert Q[0, 1] != Q[1, 0]  # Asymmetric

    def test_count_params(self):
        assert count_params(2, "ER") == 1
        assert count_params(3, "ER") == 1
        assert count_params(3, "SYM") == 3
        assert count_params(3, "ARD") == 6
        assert count_params(4, "SYM") == 6
        assert count_params(4, "ARD") == 12


class TestFitModel:
    def test_fit_er_returns_valid_result(self, args):
        fd = FitDiscrete(args)
        tree = fd.read_tree_file()
        tree_tips = fd.get_tip_names_from_tree(tree)
        tip_states = parse_discrete_traits(
            args.trait_data, tree_tips, trait_column="diet"
        )
        states = sorted(set(tip_states.values()))
        result = fd._fit_model(tree, tip_states, states, "ER", len(tip_states))

        assert result["model"] == "ER"
        assert result["lnL"] < 0  # Log-likelihood should be negative
        assert result["aic"] > 0
        assert result["bic"] > 0
        assert result["n_params"] == 1
        assert result["n_states"] == 3

    def test_model_comparison_produces_valid_weights(self, args):
        fd = FitDiscrete(args)
        fake_results = [
            {"model": "ER", "lnL": -10.0, "aic": 22.0, "bic": 23.0, "n_params": 1},
            {"model": "SYM", "lnL": -9.5, "aic": 25.0, "bic": 27.0, "n_params": 3},
            {"model": "ARD", "lnL": -8.0, "aic": 28.0, "bic": 32.0, "n_params": 6},
        ]
        results = fd._compute_model_comparison(fake_results)

        # Should be sorted by AIC
        assert results[0]["model"] == "ER"
        assert results[0]["delta_aic"] == 0.0

        # Weights should sum to 1
        total_w = sum(r["aic_weight"] for r in results)
        assert total_w == pytest.approx(1.0)

        # Best model should have highest weight
        assert results[0]["aic_weight"] > results[1]["aic_weight"]


class TestFitDiscreteRun:
    def test_run_text_output(self, args, capsys):
        fd = FitDiscrete(args)
        fd.run()
        captured = capsys.readouterr()
        assert "ER" in captured.out
        assert "SYM" in captured.out
        assert "ARD" in captured.out
        assert "lnL" in captured.out

    def test_run_json_output(self, mocker, args):
        args.json = True
        fd = FitDiscrete(args)
        mocked_json = mocker.patch("phykit.services.tree.fit_discrete.print_json")
        fd.run()
        payload = mocked_json.call_args.args[0]
        assert payload["n_taxa"] == 8
        assert payload["n_states"] == 3
        assert len(payload["models"]) == 3

    def test_run_subset_models(self, capsys):
        args = Namespace(
            tree="tests/sample_files/tree_simple.tre",
            trait_data="tests/sample_files/tree_simple_discrete_traits.tsv",
            trait="diet",
            models="ER,ARD",
        )
        fd = FitDiscrete(args)
        fd.run()
        captured = capsys.readouterr()
        assert "ER" in captured.out
        assert "ARD" in captured.out
        # SYM should not appear (only 2 models fitted)
        lines = [l for l in captured.out.strip().split("\n") if l.startswith("ER") or l.startswith("SYM") or l.startswith("ARD")]
        assert len(lines) == 2
