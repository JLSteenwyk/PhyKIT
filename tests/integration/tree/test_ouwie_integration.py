import json
import sys
from pathlib import Path

import pytest
from mock import patch

from phykit.phykit import Phykit


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_traits.tsv")
REGIMES_FILE = str(SAMPLE_FILES / "tree_simple_regimes.tsv")


@pytest.mark.integration
class TestOUwieIntegration:
    @patch("builtins.print")
    def test_basic_invocation(self, mocked_print):
        testargs = [
            "phykit",
            "ouwie",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
            "-r", REGIMES_FILE,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "OUwie Model Comparison" in all_output
        assert "BM1" in all_output

    @patch("builtins.print")
    def test_alias_fit_ouwie(self, mocked_print):
        testargs = [
            "phykit",
            "fit_ouwie",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
            "-r", REGIMES_FILE,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "OUwie Model Comparison" in all_output

    @patch("builtins.print")
    def test_alias_multi_regime_ou(self, mocked_print):
        testargs = [
            "phykit",
            "multi_regime_ou",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
            "-r", REGIMES_FILE,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "OUwie Model Comparison" in all_output

    @patch("builtins.print")
    def test_json_output(self, mocked_print):
        testargs = [
            "phykit",
            "ouwie",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
            "-r", REGIMES_FILE,
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["n_tips"] == 8
        assert "models" in payload
        assert "best_model_aicc" in payload
        assert "best_model_bic" in payload

    @patch("builtins.print")
    def test_subset_models(self, mocked_print):
        testargs = [
            "phykit",
            "ouwie",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
            "-r", REGIMES_FILE,
            "--models", "BM1,OUM,OUMVA",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert set(payload["models"].keys()) == {"BM1", "OUM", "OUMVA"}
