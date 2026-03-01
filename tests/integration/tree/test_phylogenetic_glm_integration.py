import json
import sys
from pathlib import Path

import pytest
from mock import patch

from phykit.phykit import Phykit


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
GLM_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_glm_traits.tsv")


@pytest.mark.integration
class TestPhylogeneticGLMIntegration:
    @patch("builtins.print")
    def test_binomial_basic(self, mocked_print):
        testargs = [
            "phykit",
            "phylogenetic_glm",
            "-t", TREE_SIMPLE,
            "-d", GLM_TRAITS_FILE,
            "-y", "binary_trait",
            "-x", "body_mass",
            "--family", "binomial",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Logistic MPLE" in all_output
        assert "(Intercept)" in all_output
        assert "body_mass" in all_output

    @patch("builtins.print")
    def test_poisson_basic(self, mocked_print):
        testargs = [
            "phykit",
            "phylogenetic_glm",
            "-t", TREE_SIMPLE,
            "-d", GLM_TRAITS_FILE,
            "-y", "count_trait",
            "-x", "body_mass",
            "--family", "poisson",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Poisson GEE" in all_output
        assert "(Intercept)" in all_output
        assert "Overdispersion" in all_output

    @patch("builtins.print")
    def test_alias_phylo_glm(self, mocked_print):
        testargs = [
            "phykit",
            "phylo_glm",
            "-t", TREE_SIMPLE,
            "-d", GLM_TRAITS_FILE,
            "-y", "binary_trait",
            "-x", "body_mass",
            "--family", "binomial",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Logistic MPLE" in all_output

    @patch("builtins.print")
    def test_alias_pglm(self, mocked_print):
        testargs = [
            "phykit",
            "pglm",
            "-t", TREE_SIMPLE,
            "-d", GLM_TRAITS_FILE,
            "-y", "count_trait",
            "-x", "body_mass",
            "--family", "poisson",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Poisson GEE" in all_output

    @patch("builtins.print")
    def test_json_binomial(self, mocked_print):
        testargs = [
            "phykit",
            "pglm",
            "-t", TREE_SIMPLE,
            "-d", GLM_TRAITS_FILE,
            "-y", "binary_trait",
            "-x", "body_mass",
            "--family", "binomial",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["family"] == "binomial"
        assert "alpha" in payload
        assert "coefficients" in payload
        assert payload["formula"] == "binary_trait ~ body_mass"

    @patch("builtins.print")
    def test_json_poisson(self, mocked_print):
        testargs = [
            "phykit",
            "pglm",
            "-t", TREE_SIMPLE,
            "-d", GLM_TRAITS_FILE,
            "-y", "count_trait",
            "-x", "body_mass",
            "--family", "poisson",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["family"] == "poisson"
        assert "overdispersion" in payload
        assert "coefficients" in payload

    @patch("builtins.print")
    def test_multiple_predictors(self, mocked_print):
        testargs = [
            "phykit",
            "pglm",
            "-t", TREE_SIMPLE,
            "-d", GLM_TRAITS_FILE,
            "-y", "binary_trait",
            "-x", "body_mass", "count_trait",
            "--family", "binomial",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "body_mass" in all_output
        assert "count_trait" in all_output
        assert "binary_trait ~ body_mass + count_trait" in all_output

    @patch("builtins.print")
    def test_custom_btol(self, mocked_print):
        testargs = [
            "phykit",
            "pglm",
            "-t", TREE_SIMPLE,
            "-d", GLM_TRAITS_FILE,
            "-y", "binary_trait",
            "-x", "body_mass",
            "--family", "binomial",
            "--btol", "20",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Logistic MPLE" in all_output

    @patch("builtins.print")
    def test_custom_log_alpha_bound(self, mocked_print):
        testargs = [
            "phykit",
            "pglm",
            "-t", TREE_SIMPLE,
            "-d", GLM_TRAITS_FILE,
            "-y", "binary_trait",
            "-x", "body_mass",
            "--family", "binomial",
            "--log-alpha-bound", "8",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Logistic MPLE" in all_output
