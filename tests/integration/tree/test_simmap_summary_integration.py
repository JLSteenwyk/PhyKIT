import pytest
import sys
from mock import patch
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestSimmapSummaryIntegration(object):
    @patch("builtins.print")
    def test_basic_invocation(self, mocked_print):
        testargs = [
            "phykit", "simmap_summary",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/tree_simple_discrete_traits.tsv",
            "-c", "diet", "-n", "10", "--seed", "42",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.call_count >= 1

    @patch("builtins.print")
    def test_alias_smsummary(self, mocked_print):
        testargs = [
            "phykit", "smsummary",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/tree_simple_discrete_traits.tsv",
            "-c", "diet", "-n", "10", "--seed", "42",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.call_count >= 1

    @patch("builtins.print")
    def test_alias_describe_simmap(self, mocked_print):
        testargs = [
            "phykit", "describe_simmap",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/tree_simple_discrete_traits.tsv",
            "-c", "diet", "-n", "10", "--seed", "42",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.call_count >= 1

    @patch("builtins.print")
    def test_json_output(self, mocked_print):
        testargs = [
            "phykit", "simmap_summary",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/tree_simple_discrete_traits.tsv",
            "-c", "diet", "-n", "10", "--seed", "42", "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.call_count >= 1
