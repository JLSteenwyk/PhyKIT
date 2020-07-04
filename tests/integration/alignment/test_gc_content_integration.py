import pytest
import sys
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestGCContent(object):
    @patch("builtins.print")
    def test_gc_content(self, mocked_print):
        expected_result = "0.2273"
        testargs = [
            "phykit",
            "gc_content",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]