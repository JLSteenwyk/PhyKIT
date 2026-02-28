import json
import sys
from pathlib import Path

import pytest
from mock import patch

from phykit.phykit import Phykit


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
MULTI_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_multi_traits.tsv")


@pytest.mark.integration
class TestPhylogeneticRegressionIntegration:
    @patch("builtins.print")
    def test_pgls_basic(self, mocked_print):
        testargs = [
            "phykit",
            "phylogenetic_regression",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "-y", "brain_size",
            "-x", "body_mass",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "PGLS" in all_output
        assert "(Intercept)" in all_output
        assert "body_mass" in all_output

    @patch("builtins.print")
    def test_pgls_alias(self, mocked_print):
        testargs = [
            "phykit",
            "pgls",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "-y", "brain_size",
            "-x", "body_mass",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "PGLS" in all_output

    @patch("builtins.print")
    def test_pgls_phylo_regression_alias(self, mocked_print):
        testargs = [
            "phykit",
            "phylo_regression",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "-y", "brain_size",
            "-x", "body_mass",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "PGLS" in all_output

    @patch("builtins.print")
    def test_pgls_lambda(self, mocked_print):
        testargs = [
            "phykit",
            "pgls",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "-y", "brain_size",
            "-x", "body_mass",
            "-m", "lambda",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Estimated lambda" in all_output

    @patch("builtins.print")
    def test_pgls_json(self, mocked_print):
        testargs = [
            "phykit",
            "pgls",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "-y", "brain_size",
            "-x", "body_mass",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert "coefficients" in payload
        assert "r_squared" in payload
        assert "formula" in payload
        assert payload["formula"] == "brain_size ~ body_mass"

    @patch("builtins.print")
    def test_pgls_multiple_predictors(self, mocked_print):
        testargs = [
            "phykit",
            "pgls",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "-y", "brain_size",
            "-x", "body_mass", "longevity",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "body_mass" in all_output
        assert "longevity" in all_output
        assert "brain_size ~ body_mass + longevity" in all_output

    @patch("builtins.print")
    def test_pgls_lambda_json(self, mocked_print):
        testargs = [
            "phykit",
            "pgls",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "-y", "brain_size",
            "-x", "body_mass",
            "-m", "lambda",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert "lambda" in payload
        assert "coefficients" in payload
        assert 0 <= payload["lambda"] <= 2.0
