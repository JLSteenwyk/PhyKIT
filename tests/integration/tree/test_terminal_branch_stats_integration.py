from mock import patch, call
from pathlib import Path
import pytest
import sys

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
            call(19.1996,"raccoon"),
            call(6.8004,"bear"),
            call(11.997,"sea_lion"),
            call(12.003,"seal"),
            call(100.8593,"monkey"),
            call(47.1407,"cat"),
            call(18.8795,"weasel"),
            call(25.4615,"dog"),
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
            call(0.0003,"Aspergillus_fischeri_IBT_3003"),
            call(0.0005,"Aspergillus_fischeri_IBT_3007"),
            call(0.0003,"Aspergillus_fischeri_NRRL181.GCF_000149645.1_ASM14964v1"),
            call(0.0008,"Aspergillus_fischeri_NRRL4585"),
            call(0.0017,"Aspergillus_fumigatus_Af293"),
            call(0.0005,"Aspergillus_fumigatus_CEA10"),
            call(0.0029,"Aspergillus_fumigatus_HMR_AF_270"),
            call(0.0006,"Aspergillus_fumigatus_Z5"),
            call(0.0075,"Aspergillus_oerlinghausenensis_CBS139183"),
            call(0.0004,"Aspergillus_fischeri_NRRL4161"),
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
