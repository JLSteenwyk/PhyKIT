import pytest
import sys
import json
from mock import patch, call
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestAlignmentLengthNoGaps(object):
    @patch("builtins.print")
    def test_alignment_length_invalid_command(self, mocked_print):

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 1

    @patch("builtins.print")
    def test_alignment_length_no_gaps0(self, mocked_print):
        expected_result = "3\t6\t50.0"
        testargs = [
            "phykit",
            "alignment_length_no_gaps",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_alignment_length_no_gaps1(self, mocked_print):
        expected_result = "7\t9\t77.7778"
        testargs = [
            "phykit",
            "alignment_length_no_gaps",
            f"{here.parent.parent.parent}/sample_files/test_alignment_0.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_alignment_length_no_gaps_alias0(self, mocked_print):
        expected_result = "7\t9\t77.7778"
        testargs = [
            "phykit",
            "aln_len_no_gaps",
            f"{here.parent.parent.parent}/sample_files/test_alignment_0.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_alignment_length_no_gaps_alias1(self, mocked_print):
        expected_result = "7\t9\t77.7778"
        testargs = [
            "phykit",
            "alng",
            f"{here.parent.parent.parent}/sample_files/test_alignment_0.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_alignment_length_no_gaps_incorrect_input_file(self, mocked_print):

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 1

    @patch("builtins.print")
    def test_alignment_length_no_gaps_json(self, mocked_print):
        testargs = [
            "phykit",
            "alignment_length_no_gaps",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload == {
            "alignment_length_no_gaps": 3,
            "alignment_length": 6,
            "percent_no_gaps": 50.0,
        }
