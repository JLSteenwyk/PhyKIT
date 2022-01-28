import pytest
import sys
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestSumOfPairsScore(object):
    @patch("builtins.print")
    def test_sum_of_pairs_score_full_ref(self, mocked_print):
        expected_result = 0.7714
        testargs = [
            "phykit",
            "sum_of_pairs_score",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "--reference",
            f"{here.parent.parent.parent}/sample_files/simple_reference.fa"
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_sum_of_pairs_score_short_ref(self, mocked_print):
        expected_result = 0.7714
        testargs = [
            "phykit",
            "sum_of_pairs_score",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-r",
            f"{here.parent.parent.parent}/sample_files/simple_reference.fa"
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_sum_of_pairs_score_alias0(self, mocked_print):
        expected_result = 0.7714
        testargs = [
            "phykit",
            "sops",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-r",
            f"{here.parent.parent.parent}/sample_files/simple_reference.fa"
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_sum_of_pairs_score_alias1(self, mocked_print):
        expected_result = 0.7714
        testargs = [
            "phykit",
            "sop",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-r",
            f"{here.parent.parent.parent}/sample_files/simple_reference.fa"
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]