import json
import os
import sys
import tempfile
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
class TestRateHeterogeneityIntegration:
    @patch("builtins.print")
    def test_basic_invocation(self, mocked_print):
        testargs = [
            "phykit",
            "rate_heterogeneity",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
            "-r", REGIMES_FILE,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Rate Heterogeneity Test" in all_output
        assert "Single-rate model" in all_output
        assert "Multi-rate model" in all_output

    @patch("builtins.print")
    def test_alias_brownie(self, mocked_print):
        testargs = [
            "phykit",
            "brownie",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
            "-r", REGIMES_FILE,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Rate Heterogeneity Test" in all_output

    @patch("builtins.print")
    def test_alias_rh(self, mocked_print):
        testargs = [
            "phykit",
            "rh",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
            "-r", REGIMES_FILE,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Rate Heterogeneity Test" in all_output

    @patch("builtins.print")
    def test_json_output(self, mocked_print):
        testargs = [
            "phykit",
            "rh",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
            "-r", REGIMES_FILE,
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["n_tips"] == 8
        assert "single_rate" in payload
        assert "multi_rate" in payload
        assert "lrt" in payload

    @patch("builtins.print")
    def test_with_nsim_and_seed(self, mocked_print):
        testargs = [
            "phykit",
            "rh",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
            "-r", REGIMES_FILE,
            "-n", "10",
            "--seed", "42",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["lrt"]["n_sim"] == 10
        assert payload["lrt"]["sim_p_value"] is not None

    def test_plot_output(self):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            testargs = [
                "phykit",
                "rh",
                "-t", TREE_SIMPLE,
                "-d", TRAITS_FILE,
                "-r", REGIMES_FILE,
                "--plot", tmppath,
            ]
            with patch.object(sys, "argv", testargs):
                Phykit()
            assert os.path.exists(tmppath)
            assert os.path.getsize(tmppath) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)
