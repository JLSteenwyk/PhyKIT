import pytest
import sys
from math import isclose
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestTotalTreeLength(object):
    @patch("builtins.print")
    def test_total_tree_length0(self, mocked_print):
        expected_result = 277.2772
        testargs = [
            "phykit",
            "total_tree_length",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_total_tree_length1(self, mocked_print):
        expected_result = 0.0675
        testargs = [
            "phykit",
            "tree_len",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_total_tree_length_alias(self, mocked_print):
        expected_result = 277.2772
        testargs = [
            "phykit",
            "tree_len",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_total_tree_length_incorrect_file_path(self, mocked_print):
        testargs = [
            "phykit",
            "total_tree_length",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tr",
        ]

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2
