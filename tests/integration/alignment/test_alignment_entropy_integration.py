import pytest
import sys
import json
from mock import patch, call
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestAlignmentEntropy(object):
    @patch("builtins.print")
    def test_alignment_entropy0(self, mocked_print):
        expected_result = 0.657
        testargs = [
            "phykit",
            "alignment_entropy",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_alignment_entropy_alias(self, mocked_print):
        expected_result = 0.657
        testargs = [
            "phykit",
            "aln_entropy",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_alignment_entropy_incorrect_input_file(self, mocked_print):
        testargs = [
            "phykit",
            "alignment_entropy",
            f"{here.parent.parent.parent}/sample_files/test_alignment_1.f",
        ]
        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_alignment_entropy_json_summary(self, mocked_print):
        testargs = [
            "phykit",
            "alignment_entropy",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload == {"mean_entropy": 0.657, "verbose": False}

    @patch("builtins.print")
    def test_alignment_entropy_json_verbose(self, mocked_print):
        testargs = [
            "phykit",
            "alignment_entropy",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-v",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["verbose"] is True
        assert payload["rows"][0] == payload["sites"][0]
        assert payload["sites"][0] == {"site": 1, "entropy": 0.0}
