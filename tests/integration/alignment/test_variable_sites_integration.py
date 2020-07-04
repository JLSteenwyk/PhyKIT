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
    def test_variable_sites(self, mocked_print):
        expected_result = "4\t6\t66.6667"
        testargs = [
            "phykit",
            "variable_sites",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]
