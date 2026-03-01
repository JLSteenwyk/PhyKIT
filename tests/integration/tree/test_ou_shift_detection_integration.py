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


@pytest.mark.integration
class TestOUShiftDetectionIntegration:
    @patch("builtins.print")
    def test_basic_invocation(self, mocked_print):
        testargs = [
            "phykit",
            "ou_shift_detection",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "OU Shift Detection" in all_output

    @patch("builtins.print")
    def test_alias_l1ou(self, mocked_print):
        testargs = [
            "phykit",
            "l1ou",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "OU Shift Detection" in all_output

    @patch("builtins.print")
    def test_alias_ou_shifts(self, mocked_print):
        testargs = [
            "phykit",
            "ou_shifts",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "OU Shift Detection" in all_output

    @patch("builtins.print")
    def test_alias_detect_shifts(self, mocked_print):
        testargs = [
            "phykit",
            "detect_shifts",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "OU Shift Detection" in all_output

    @patch("builtins.print")
    def test_json_output(self, mocked_print):
        testargs = [
            "phykit",
            "l1ou",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["n_tips"] == 8
        assert "n_shifts" in payload
        assert "alpha" in payload
        assert "sigma2" in payload
        assert "theta_root" in payload
        assert "log_likelihood" in payload
        assert "shifts" in payload

    @patch("builtins.print")
    def test_criterion_bic(self, mocked_print):
        testargs = [
            "phykit",
            "l1ou",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
            "--criterion", "BIC",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "OU Shift Detection" in all_output

    @patch("builtins.print")
    def test_criterion_aicc(self, mocked_print):
        testargs = [
            "phykit",
            "l1ou",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
            "--criterion", "AICc",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "OU Shift Detection" in all_output

    @patch("builtins.print")
    def test_max_shifts(self, mocked_print):
        testargs = [
            "phykit",
            "l1ou",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
            "--max-shifts", "2",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "OU Shift Detection" in all_output

    @patch("builtins.print")
    def test_json_with_max_shifts(self, mocked_print):
        testargs = [
            "phykit",
            "l1ou",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
            "--max-shifts", "1",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["n_shifts"] <= 1
