import pytest
import sys
from math import isclose
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestIBS(object):
    @patch("builtins.print")
    def test_internal_branch_stats0(self, mocked_print):
        testargs = [
            "phykit",
            "internal_branch_stats",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call("mean: 6.9872"),
            call("median: 3.8738"),
            call("25th percentile: 2.0946"),
            call("75th percentile: 7.5297"),
            call("minimum: 0.846"),
            call("maximum: 20.592"),
            call("standard deviation: 8.0114"),
            call("variance: 64.1826")
        ]

    @patch("builtins.print")
    def test_internal_branch_stats1(self, mocked_print):
        testargs = [
            "phykit",
            "internal_branch_stats",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call("mean: 0.0074"),
            call("median: 0.0007"),
            call("25th percentile: 0.0004"),
            call("75th percentile: 0.0083"),
            call("minimum: 0.0002"),
            call("maximum: 0.0337"),
            call("standard deviation: 0.0129"),
            call("variance: 0.0002")
        ]

    @patch("builtins.print")
    def test_internal_branch_stats_verbose0(self, mocked_print):
        testargs = [
            "phykit",
            "internal_branch_stats",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-v"
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call(0.846),
            call(3.8738),
            call(7.5297),
            call(2.0946),
            call(20.592)
        ]

    @patch("builtins.print")
    def test_internal_branch_verbose1(self, mocked_print):
        testargs = [
            "phykit",
            "internal_branch_stats",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            "-v"
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call(0.0002),
            call(0.0002),
            call(0.0006),
            call(0.0158),
            call(0.0337),
            call(0.0008),
            call(0.0007)
        ]

    @patch("builtins.print")
    def test_internal_branch_stats_alias(self, mocked_print):
        testargs = [
            "phykit",
            "ibs",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call("mean: 6.9872"),
            call("median: 3.8738"),
            call("25th percentile: 2.0946"),
            call("75th percentile: 7.5297"),
            call("minimum: 0.846"),
            call("maximum: 20.592"),
            call("standard deviation: 8.0114"),
            call("variance: 64.1826")
        ]

    @patch("builtins.print")
    def test_internal_branch_stats_alias(self, mocked_print):
        testargs = [
            "phykit",
            "internal_branch_stats",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tr",
        ]
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2