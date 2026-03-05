import pytest
import sys
import json
from mock import patch, call
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestVariableSites(object):
    @patch("builtins.print")
    def test_variable_sites0(self, mocked_print):
        expected_result = "4\t6\t66.6667"
        testargs = [
            "phykit",
            "variable_sites",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_variable_sites1(self, mocked_print):
        expected_result = "2\t9\t22.2222"
        testargs = [
            "phykit",
            "variable_sites",
            f"{here.parent.parent.parent}/sample_files/test_alignment_0.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_variable_sites2(self, mocked_print):
        expected_result = "1\t3\t33.3333"
        testargs = [
            "phykit",
            "variable_sites",
            f"{here.parent.parent.parent}/sample_files/test_alignment_1.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_variable_sites_alias(self, mocked_print):
        expected_result = "1\t3\t33.3333"
        testargs = [
            "phykit",
            "vs",
            f"{here.parent.parent.parent}/sample_files/test_alignment_1.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_variable_sites_incorrect_input_file(self, mocked_print):
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_variable_sites_json(self, mocked_print):
        testargs = [
            "phykit",
            "variable_sites",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload == {
            "variable_sites": 4,
            "alignment_length": 6,
            "percent_variable_sites": 66.6667,
        }
