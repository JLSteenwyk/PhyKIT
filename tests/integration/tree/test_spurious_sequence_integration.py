from mock import patch, call
from pathlib import Path
import pytest
import sys
import json

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestSpuriousSequence(object):
    @patch("builtins.print")
    def test_spurious_sequence(self, mocked_print):
        expected_result = "None"
        testargs = [
            "phykit",
            "spurious_sequence",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_spurious_sequence_alias0(self, mocked_print):
        expected_result = "None"
        testargs = [
            "phykit",
            "spurious_seq",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_spurious_sequence_alias1(self, mocked_print):
        expected_result = "None"
        testargs = [
            "phykit",
            "ss",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_spurious_sequence_custom_factor(self, mocked_print):
        # With factor 2 and terminal branch median of ~19.04, both monkey and cat exceed threshold
        expected_results = [
            call("monkey\t100.8593\t38.0791\t19.0396"),
            call("cat\t47.1407\t38.0791\t19.0396")
        ]
        testargs = [
            "phykit",
            "spurious_sequence",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-f",
            "2"
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == expected_results

    @patch("builtins.print")
    def test_spurious_sequence_incorrect_file_path(self, mocked_print):
        testargs = [
            "phykit",
            "spurious_sequence",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tr",  # Invalid file path
        ]
        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_spurious_sequence_json_none(self, mocked_print):
        testargs = [
            "phykit",
            "spurious_sequence",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload == {"rows": [], "spurious_sequences": []}

    @patch("builtins.print")
    def test_spurious_sequence_json_custom_factor(self, mocked_print):
        testargs = [
            "phykit",
            "spurious_sequence",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-f",
            "2",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["rows"][0] == payload["spurious_sequences"][0]
        assert payload["spurious_sequences"][0] == {
            "taxon": "monkey",
            "branch_length": 100.8593,
            "threshold": 38.0791,
            "median": 19.0396,
        }
