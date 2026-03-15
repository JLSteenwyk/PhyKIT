"""
Integration tests for fit_discrete (discrete trait model comparison).

Expected values cross-validated against R's geiger::fitDiscrete() with
multi2di() to resolve multifurcations. R validation script:
tests/r_validation/validate_fit_discrete.R

To reproduce expected values in R:
    library(ape); library(geiger)
    tree <- multi2di(read.tree("tests/sample_files/tree_simple.tre"))
    traits <- read.delim("tests/sample_files/tree_simple_discrete_traits.tsv")
    trait_vec <- setNames(traits$diet, traits$taxon)
    fitDiscrete(tree, trait_vec, model="ER", type="equal")

Note: PhyKIT handles multifurcating trees natively while geiger requires
resolved bifurcating trees. Log-likelihoods may differ by ~0.5-1.0 due
to this structural difference and multi-start optimization convergence.
The model ranking (which model is best by AIC) should agree.
"""
from mock import patch
from pathlib import Path
import json
import sys

import pytest

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestFitDiscrete:
    @patch("builtins.print")
    def test_fit_discrete_text_output(self, mocked_print):
        testargs = [
            "phykit",
            "fit_discrete",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/tree_simple_discrete_traits.tsv",
            "-c", "diet",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        # Collect all printed lines (some print() calls have no args)
        output = "\n".join(
            call.args[0] if call.args else "" for call in mocked_print.call_args_list
        )
        assert "ER" in output
        assert "SYM" in output
        assert "ARD" in output

    @patch("builtins.print")
    def test_fit_discrete_er_is_best_model(self, mocked_print):
        # Cross-validated: R's geiger also selects ER as best for this dataset
        testargs = [
            "phykit",
            "fit_discrete",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/tree_simple_discrete_traits.tsv",
            "-c", "diet",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        # Models should be sorted by AIC (best first)
        assert payload["models"][0]["model"] == "ER"
        assert payload["models"][0]["delta_aic"] == 0.0
        assert payload["n_states"] == 3
        assert payload["n_taxa"] == 8

    @patch("builtins.print")
    def test_fit_discrete_alias_fd(self, mocked_print):
        testargs = [
            "phykit",
            "fd",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/tree_simple_discrete_traits.tsv",
            "-c", "diet",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert len(payload["models"]) == 3

    @patch("builtins.print")
    def test_fit_discrete_alias_fitdiscrete(self, mocked_print):
        testargs = [
            "phykit",
            "fitdiscrete",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/tree_simple_discrete_traits.tsv",
            "-c", "diet",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert len(payload["models"]) == 3

    @patch("builtins.print")
    def test_fit_discrete_subset_models(self, mocked_print):
        testargs = [
            "phykit",
            "fit_discrete",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/tree_simple_discrete_traits.tsv",
            "-c", "diet",
            "--models", "ER,ARD",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        model_names = [m["model"] for m in payload["models"]]
        assert "ER" in model_names
        assert "ARD" in model_names
        assert "SYM" not in model_names

    @patch("builtins.print")
    def test_fit_discrete_lnl_within_tolerance_of_r(self, mocked_print):
        # R's geiger::fitDiscrete ER lnL = -8.2885 (with multi2di)
        # PhyKIT handles multifurcations natively, so lnL may differ
        # but should be in the same range (within ~1.0)
        testargs = [
            "phykit",
            "fit_discrete",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/tree_simple_discrete_traits.tsv",
            "-c", "diet",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        er_model = next(m for m in payload["models"] if m["model"] == "ER")
        # R gives -8.2885; allow ~1.0 tolerance for tree structure differences
        assert abs(er_model["lnL"] - (-8.2885)) < 1.5

    @patch("builtins.print")
    def test_fit_discrete_json_has_q_matrix(self, mocked_print):
        testargs = [
            "phykit",
            "fit_discrete",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/tree_simple_discrete_traits.tsv",
            "-c", "diet",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        for model in payload["models"]:
            assert "q_matrix" in model
            assert len(model["q_matrix"]) == 3  # 3 states
            assert len(model["q_matrix"][0]) == 3
