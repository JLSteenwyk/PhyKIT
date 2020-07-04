import pytest
import sys
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestPairwiseIdentity(object):
    @patch("builtins.print")
    def test_pairwise_identity(self, mocked_print):
        testargs = [
            "phykit",
            "pairwise_identity",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call("mean: 0.4667"),
            call("median: 0.5"),
            call("25th percentile: 0.3333"),
            call("75th percentile: 0.625"),
            call("minimum: 0.1667"),
            call("maximum: 0.8333"),
            call("standard deviation: 0.2194"),
            call("variance: 0.0481")
        ]