import pytest
import sys
from mock import patch, call
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestColumnScore(object):
    @patch("builtins.print")
    def test_column_score_full_ref(self, mocked_print):
        expected_result = 0.8333
        testargs = [
            "phykit",
            "column_score",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "--reference",
            f"{here.parent.parent.parent}/sample_files/simple_reference.fa"
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_column_score_short_ref(self, mocked_print):
        expected_result = 0.8333
        testargs = [
            "phykit",
            "column_score",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-r",
            f"{here.parent.parent.parent}/sample_files/simple_reference.fa"
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_column_score_alias(self, mocked_print):
        expected_result = 0.8333
        testargs = [
            "phykit",
            "cs",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-r",
            f"{here.parent.parent.parent}/sample_files/simple_reference.fa"
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]
