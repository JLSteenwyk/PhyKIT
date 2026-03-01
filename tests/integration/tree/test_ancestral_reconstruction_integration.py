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
MULTI_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_multi_traits.tsv")


@pytest.mark.integration
class TestAncestralReconstructionIntegration:
    @patch("builtins.print")
    def test_basic_fast(self, mocked_print):
        testargs = [
            "phykit",
            "ancestral_state_reconstruction",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Ancestral State Reconstruction" in all_output
        assert "fast" in all_output
        assert "N1 (root)" in all_output

    @patch("builtins.print")
    def test_alias_asr(self, mocked_print):
        testargs = [
            "phykit",
            "asr",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Ancestral State Reconstruction" in all_output

    @patch("builtins.print")
    def test_alias_anc_recon(self, mocked_print):
        testargs = [
            "phykit",
            "anc_recon",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Ancestral State Reconstruction" in all_output

    @patch("builtins.print")
    def test_ml_method(self, mocked_print):
        testargs = [
            "phykit",
            "asr",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
            "-m", "ml",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "ml" in all_output

    @patch("builtins.print")
    def test_ci_flag(self, mocked_print):
        testargs = [
            "phykit",
            "asr",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
            "--ci",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "95% CI" in all_output

    @patch("builtins.print")
    def test_json_output(self, mocked_print):
        testargs = [
            "phykit",
            "asr",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert "method" in payload
        assert payload["method"] == "fast"
        assert "ancestral_estimates" in payload
        assert "tip_values" in payload

    @patch("builtins.print")
    def test_plot_output(self, mocked_print):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            plot_path = f.name

        try:
            testargs = [
                "phykit",
                "asr",
                "-t", TREE_SIMPLE,
                "-d", TRAITS_FILE,
                "--plot", plot_path,
            ]
            with patch.object(sys, "argv", testargs):
                Phykit()

            assert os.path.exists(plot_path)
            assert os.path.getsize(plot_path) > 0
        finally:
            if os.path.exists(plot_path):
                os.unlink(plot_path)

    @patch("builtins.print")
    def test_multi_trait_with_c_flag(self, mocked_print):
        testargs = [
            "phykit",
            "asr",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "-c", "body_mass",
            "--ci",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "body_mass" in all_output
        assert "95% CI" in all_output
