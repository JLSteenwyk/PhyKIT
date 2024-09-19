from mock import patch, call
from pathlib import Path
import pytest
import sys

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
        expected_result = "monkey\t100.8593\t50.9231\t25.4615"
        testargs = [
            "phykit",
            "spurious_sequence",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-f",
            "2"
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_spurious_sequence_incorrect_file_path(self, mocked_print):
        testargs = [
            "phykit",
            "spurious_sequence",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
        ]
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2
