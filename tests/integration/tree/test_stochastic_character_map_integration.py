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
DISCRETE_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_discrete_traits.tsv")


@pytest.mark.integration
class TestStochasticCharacterMapIntegration:
    @patch("builtins.print")
    def test_basic_er(self, mocked_print):
        testargs = [
            "phykit",
            "stochastic_character_map",
            "-t", TREE_SIMPLE,
            "-d", DISCRETE_TRAITS_FILE,
            "-c", "diet",
            "-n", "10",
            "--seed", "42",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "SIMMAP" in all_output
        assert "Model: ER" in all_output

    @patch("builtins.print")
    def test_alias_simmap(self, mocked_print):
        testargs = [
            "phykit",
            "simmap",
            "-t", TREE_SIMPLE,
            "-d", DISCRETE_TRAITS_FILE,
            "-c", "diet",
            "-n", "10",
            "--seed", "42",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "SIMMAP" in all_output

    @patch("builtins.print")
    def test_alias_scm(self, mocked_print):
        testargs = [
            "phykit",
            "scm",
            "-t", TREE_SIMPLE,
            "-d", DISCRETE_TRAITS_FILE,
            "-c", "diet",
            "-n", "10",
            "--seed", "42",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "SIMMAP" in all_output

    @patch("builtins.print")
    def test_sym_model(self, mocked_print):
        testargs = [
            "phykit",
            "simmap",
            "-t", TREE_SIMPLE,
            "-d", DISCRETE_TRAITS_FILE,
            "-c", "diet",
            "-m", "SYM",
            "-n", "10",
            "--seed", "42",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Model: SYM" in all_output

    @patch("builtins.print")
    def test_ard_model(self, mocked_print):
        testargs = [
            "phykit",
            "simmap",
            "-t", TREE_SIMPLE,
            "-d", DISCRETE_TRAITS_FILE,
            "-c", "diet",
            "-m", "ARD",
            "-n", "10",
            "--seed", "42",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Model: ARD" in all_output

    @patch("builtins.print")
    def test_json_output(self, mocked_print):
        testargs = [
            "phykit",
            "simmap",
            "-t", TREE_SIMPLE,
            "-d", DISCRETE_TRAITS_FILE,
            "-c", "diet",
            "-n", "10",
            "--seed", "42",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["model"] == "ER"
        assert payload["nsim"] == 10
        assert "q_matrix" in payload
        assert "log_likelihood" in payload
        assert "states" in payload
        assert "mean_dwelling_times" in payload

    @patch("builtins.print")
    def test_reproducibility_with_seed(self, mocked_print):
        testargs = [
            "phykit",
            "simmap",
            "-t", TREE_SIMPLE,
            "-d", DISCRETE_TRAITS_FILE,
            "-c", "diet",
            "-n", "10",
            "--seed", "123",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload1 = json.loads(mocked_print.call_args.args[0])

        mocked_print.reset_mock()

        with patch.object(sys, "argv", testargs):
            Phykit()
        payload2 = json.loads(mocked_print.call_args.args[0])

        # Same seed should give same results
        assert payload1["log_likelihood"] == payload2["log_likelihood"]
        assert payload1["mean_dwelling_times"] == payload2["mean_dwelling_times"]

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
                "simmap",
                "-t", TREE_SIMPLE,
                "-d", DISCRETE_TRAITS_FILE,
                "-c", "diet",
                "-n", "10",
                "--seed", "42",
                "--plot", tmppath,
            ]
            with patch.object(sys, "argv", testargs):
                Phykit()
            assert os.path.exists(tmppath)
            assert os.path.getsize(tmppath) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    def test_stochastic_character_map_circular(self):
        """--circular flag produces a circular layout stochastic map plot."""
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            testargs = [
                "phykit",
                "simmap",
                "-t", TREE_SIMPLE,
                "-d", DISCRETE_TRAITS_FILE,
                "-c", "diet",
                "-n", "10",
                "--seed", "42",
                "--plot", tmppath,
                "--circular",
            ]
            with patch.object(sys, "argv", testargs):
                Phykit()
            assert os.path.exists(tmppath)
            assert os.path.getsize(tmppath) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)
