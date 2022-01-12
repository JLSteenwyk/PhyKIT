import pytest
import sys
from math import isclose
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestTBS(object):
    @patch("builtins.print")
    def test_terminal_branch_stats0(self, mocked_print):
        testargs = [
            "phykit",
            "terminal_branch_stats",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call("mean: 30.2926"),
            call("median: 19.0396"),
            call("25th percentile: 12.0015"),
            call("75th percentile: 30.8813"),
            call("minimum: 6.8004"),
            call("maximum: 100.8593"),
            call("standard deviation: 31.0789"),
            call("variance: 965.8987")
        ]

    @patch("builtins.print")
    def test_terminal_branch_stats1(self, mocked_print):
        testargs = [
            "phykit",
            "terminal_branch_stats",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call("mean: 0.0016"),
            call("median: 0.0005"),
            call("25th percentile: 0.0004"),
            call("75th percentile: 0.0015"),
            call("minimum: 0.0003"),
            call("maximum: 0.0075"),
            call("standard deviation: 0.0023"),
            call("variance: 0.0")
        ]

    @patch("builtins.print")
    def test_terminal_branch_stats_verbose0(self, mocked_print):
        testargs = [
            "phykit",
            "terminal_branch_stats",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-v"
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call(19.1996),
            call(6.8004),
            call(11.997),
            call(12.003),
            call(100.8593),
            call(47.1407),
            call(18.8795),
            call(25.4615),
        ]

    @patch("builtins.print")
    def test_terminal_branch_verbose1(self, mocked_print):
        testargs = [
            "phykit",
            "terminal_branch_stats",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            "-v"
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call(0.0003),
            call(0.0005),
            call(0.0003),
            call(0.0008),
            call(0.0017),
            call(0.0005),
            call(0.0029),
            call(0.0006),
            call(0.0075),
            call(0.0004),
        ]

    @patch("builtins.print")
    def test_terminal_branch_stats_alias(self, mocked_print):
        testargs = [
            "phykit",
            "tbs",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call("mean: 30.2926"),
            call("median: 19.0396"),
            call("25th percentile: 12.0015"),
            call("75th percentile: 30.8813"),
            call("minimum: 6.8004"),
            call("maximum: 100.8593"),
            call("standard deviation: 31.0789"),
            call("variance: 965.8987")
        ]

    @patch("builtins.print")
    def test_terminal_branch_stats_alias(self, mocked_print):
        testargs = [
            "phykit",
            "terminal_branch_stats",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tr",
        ]
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2