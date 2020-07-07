import pytest
import sys
from math import isclose
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestTree(object):
    @patch("builtins.print")
    def test_treeness_over_rcv(self, mocked_print):
        expected_result = "0.35\t0.126\t0.36"
        testargs = [
            "phykit",
            "treeness_over_rcv",
            "-t",
            f"{here.parent.parent}/sample_files/tree_simple.tre",
            "-a",
            f"{here.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]