import pytest
import sys
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestAlignmentLength(object):
    @patch("builtins.print")
    def test_alignment_length_incorrect_file_path(self, mocked_print):
        expected_result = "Input file could not be read. Please check input file argument."
        testargs = [
            "phykit",
            "alignment_length",
            f"whoa",
        ]

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_alignment_length_test0(self, mocked_print):
        expected_result = 6
        testargs = [
            "phykit",
            "alignment_length",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_alignment_length_test1(self, mocked_print):
        expected_result = 9
        testargs = [
            "phykit",
            "alignment_length",
            f"{here.parent.parent.parent}/sample_files/test_alignment_0.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_alignment_length_test_alias0(self, mocked_print):
        expected_result = 9
        testargs = [
            "phykit",
            "aln_len",
            f"{here.parent.parent.parent}/sample_files/test_alignment_0.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_alignment_length_test_alias1(self, mocked_print):
        expected_result = 9
        testargs = [
            "phykit",
            "al",
            f"{here.parent.parent.parent}/sample_files/test_alignment_0.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_alignment_length_incorrect_input_file(self, mocked_print):
        expected_result = "Input file could not be read. Please check input file argument."
        testargs = [
            "phykit",
            "alignment_length",
            f"{here.parent.parent.parent}/sample_files/test_trees.txt",
        ]

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2