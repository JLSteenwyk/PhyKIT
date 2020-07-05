import pytest
import sys
from mock import patch, call
from pathlib import Path
from textwrap import dedent

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
        testargs = [
            "phykit",
            "rcv",
            f"{here.parent.parent.parent}/sample_files/test_alignment_1.f",
        ]
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2