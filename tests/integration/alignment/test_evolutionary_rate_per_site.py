import pytest
import sys
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestEvolutionaryRatePerSite(object):
    @patch("builtins.print")
    def test_evolutionary_rate_per_site(self, mocked_print):
        expected_result_0 = dedent(
            """1\t0.0"""
        )
        expected_result_1 = dedent(
            """2\t0.5"""
        )
        expected_result_2 = dedent(
            """3\t0.48"""  
        )
        expected_result_3 = dedent(
            """4\t0.0"""
        )
        expected_result_4 = dedent(
            """5\t0.48"""
        )
        expected_result_5 = dedent(
            """6\t0.5"""
        )
        testargs = [
            "phykit",
            "evolutionary_rate_per_site",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
            call(expected_result_2),
            call(expected_result_3),
            call(expected_result_4),
            call(expected_result_5),
        ]

    @patch("builtins.print")
    def test_evolutionary_rate_per_site_alias0(self, mocked_print):
        expected_result_0 = dedent(
            """1\t0.0"""
        )
        expected_result_1 = dedent(
            """2\t0.5"""
        )
        expected_result_2 = dedent(
            """3\t0.48"""  
        )
        expected_result_3 = dedent(
            """4\t0.0"""
        )
        expected_result_4 = dedent(
            """5\t0.48"""
        )
        expected_result_5 = dedent(
            """6\t0.5"""
        )
        testargs = [
            "phykit",
            "evo_rate_per_site",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
            call(expected_result_2),
            call(expected_result_3),
            call(expected_result_4),
            call(expected_result_5),
        ]

    @patch("builtins.print")
    def test_evolutionary_rate_per_site_alias1(self, mocked_print):
        expected_result_0 = dedent(
            """1\t0.0"""
        )
        expected_result_1 = dedent(
            """2\t0.5"""
        )
        expected_result_2 = dedent(
            """3\t0.48"""  
        )
        expected_result_3 = dedent(
            """4\t0.0"""
        )
        expected_result_4 = dedent(
            """5\t0.48"""
        )
        expected_result_5 = dedent(
            """6\t0.5"""
        )
        testargs = [
            "phykit",
            "erps",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
            call(expected_result_2),
            call(expected_result_3),
            call(expected_result_4),
            call(expected_result_5),
        ]
