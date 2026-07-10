from mock import patch, call
from pathlib import Path
import pytest
import sys
import json

from phykit.phykit import Phykit

here = Path(__file__)


def summary_call(
    mean,
    median,
    twenty_fifth,
    seventy_fifth,
    minimum,
    maximum,
    standard_deviation,
    variance,
):
    return call(
        (
            f"mean: {mean}\n"
            f"median: {median}\n"
            f"25th percentile: {twenty_fifth}\n"
            f"75th percentile: {seventy_fifth}\n"
            f"minimum: {minimum}\n"
            f"maximum: {maximum}\n"
            f"standard deviation: {standard_deviation}\n"
            f"variance: {variance}"
        )
    )


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
            summary_call(
                "30.2926", "19.0396", "12.0015", "30.8813",
                "6.8004", "100.8593", "31.0789", "965.8987",
            )
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
            summary_call(
                "0.0016", "0.0005", "0.0004", "0.0015",
                "0.0003", "0.0075", "0.0023", "0.0",
            )
        ]

    @patch("builtins.print")
    def test_terminal_branch_stats_verbose0(self, mocked_print):
        expected_result = (
            "19.1996 raccoon\n"
            "6.8004 bear\n"
            "11.997 sea_lion\n"
            "12.003 seal\n"
            "100.8593 monkey\n"
            "47.1407 cat\n"
            "18.8795 weasel\n"
            "25.4615 dog"
        )
        testargs = [
            "phykit",
            "terminal_branch_stats",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-v"
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_terminal_branch_verbose1(self, mocked_print):
        expected_result = (
            "0.0003 Aspergillus_fischeri_IBT_3003\n"
            "0.0005 Aspergillus_fischeri_IBT_3007\n"
            "0.0003 Aspergillus_fischeri_NRRL181.GCF_000149645.1_ASM14964v1\n"
            "0.0008 Aspergillus_fischeri_NRRL4585\n"
            "0.0017 Aspergillus_fumigatus_Af293\n"
            "0.0005 Aspergillus_fumigatus_CEA10\n"
            "0.0029 Aspergillus_fumigatus_HMR_AF_270\n"
            "0.0006 Aspergillus_fumigatus_Z5\n"
            "0.0075 Aspergillus_oerlinghausenensis_CBS139183\n"
            "0.0004 Aspergillus_fischeri_NRRL4161"
        )
        testargs = [
            "phykit",
            "terminal_branch_stats",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            "-v"
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_terminal_branch_stats_missing_args(self, mocked_print):
        testargs = [
            "phykit",
            "tbs",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            summary_call(
                "30.2926", "19.0396", "12.0015", "30.8813",
                "6.8004", "100.8593", "31.0789", "965.8987",
            )
        ]

    @patch("builtins.print")
    def test_terminal_branch_stats_alias(self, mocked_print):
        testargs = ["phykit", "tbs"]
        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_terminal_branch_stats_json_summary(self, mocked_print):
        testargs = [
            "phykit",
            "terminal_branch_stats",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["verbose"] is False
        assert round(payload["summary"]["mean"], 4) == 30.2926

    @patch("builtins.print")
    def test_terminal_branch_stats_json_verbose(self, mocked_print):
        testargs = [
            "phykit",
            "terminal_branch_stats",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-v",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["verbose"] is True
        assert payload["rows"][0] == payload["tips"][0]
        assert payload["tips"][0] == {"length": 19.1996, "taxon": "raccoon"}
