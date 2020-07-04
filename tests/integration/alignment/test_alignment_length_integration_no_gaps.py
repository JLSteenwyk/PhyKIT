import pytest
import sys
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestAlignmentLengthNoGaps(object):
    @patch("builtins.print")
    def test_alignment_length_no_gaps(self, mocked_print):
        expected_result = "3\t6\t50.0"
        testargs = [
            "phykit",
            "alignment_length_no_gaps",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]