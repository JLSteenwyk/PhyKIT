import pytest
import sys
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
