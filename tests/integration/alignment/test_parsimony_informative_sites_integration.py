import pytest
import sys
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestParsimonyInformativeSites(object):
    @patch("builtins.print")
    def test_parsimony_informative_sites(self, mocked_print):
        expected_result = "3\t6\t50.0"
        testargs = [
            "phykit",
            "parsimony_informative_sites",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]