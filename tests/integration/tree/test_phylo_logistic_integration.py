import json
import sys
from pathlib import Path

import pytest
from mock import patch

from phykit.phykit import Phykit


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
BINARY_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_binary_traits.tsv")
GLM_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_glm_traits.tsv")


@pytest.mark.integration
class TestPhyloLogisticIntegration:
    @patch("builtins.print")
    def test_basic_invocation(self, mocked_print):
        testargs = [
            "phykit",
            "phylo_logistic",
            "-t", TREE_SIMPLE,
            "-d", GLM_TRAITS_FILE,
            "--response", "binary_trait",
            "--predictor", "body_mass",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Phylogenetic Logistic Regression" in all_output
        assert "(Intercept)" in all_output
        assert "body_mass" in all_output

    @patch("builtins.print")
    def test_alias_phylo_logreg(self, mocked_print):
        testargs = [
            "phykit",
            "phylo_logreg",
            "-t", TREE_SIMPLE,
            "-d", GLM_TRAITS_FILE,
            "--response", "binary_trait",
            "--predictor", "body_mass",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Phylogenetic Logistic Regression" in all_output

    @patch("builtins.print")
    def test_alias_plogreg(self, mocked_print):
        testargs = [
            "phykit",
            "plogreg",
            "-t", TREE_SIMPLE,
            "-d", GLM_TRAITS_FILE,
            "--response", "binary_trait",
            "--predictor", "body_mass",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Phylogenetic Logistic Regression" in all_output

    @patch("builtins.print")
    def test_json_output(self, mocked_print):
        testargs = [
            "phykit",
            "plogreg",
            "-t", TREE_SIMPLE,
            "-d", GLM_TRAITS_FILE,
            "--response", "binary_trait",
            "--predictor", "body_mass",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["method"] == "logistic_MPLE"
        assert "alpha" in payload
        assert "coefficients" in payload
        assert "(Intercept)" in payload["coefficients"]
        assert payload["response"] == "binary_trait"

    @patch("builtins.print")
    def test_with_binary_traits_file(self, mocked_print):
        testargs = [
            "phykit",
            "phylo_logistic",
            "-t", TREE_SIMPLE,
            "-d", BINARY_TRAITS_FILE,
            "--response", "has_wings",
            "--predictor", "body_mass",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Phylogenetic Logistic Regression" in all_output
        assert "has_wings" in all_output or "body_mass" in all_output

    @patch("builtins.print")
    def test_multiple_predictors_comma_separated(self, mocked_print):
        testargs = [
            "phykit",
            "phylo_logistic",
            "-t", TREE_SIMPLE,
            "-d", GLM_TRAITS_FILE,
            "--response", "binary_trait",
            "--predictor", "body_mass,count_trait",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "body_mass" in all_output
        assert "count_trait" in all_output
