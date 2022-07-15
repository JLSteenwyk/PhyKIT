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
            call(0.846, "raccoon;bear"), 
            call(3.8738, "sea_lion;seal;monkey;cat;weasel"),
            call(7.5297, "sea_lion;seal"),
            call(2.0946, "monkey;cat;weasel"),
            call(20.592, "monkey;cat")
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
            call(0.0002,"Aspergillus_fischeri_IBT_3007;Aspergillus_fischeri_NRRL181.GCF_000149645.1_ASM14964v1;Aspergillus_fischeri_NRRL4585;Aspergillus_fumigatus_Af293;Aspergillus_fumigatus_CEA10;Aspergillus_fumigatus_HMR_AF_270;Aspergillus_fumigatus_Z5;Aspergillus_oerlinghausenensis_CBS139183"),
            call(0.0002,"Aspergillus_fischeri_NRRL181.GCF_000149645.1_ASM14964v1;Aspergillus_fischeri_NRRL4585;Aspergillus_fumigatus_Af293;Aspergillus_fumigatus_CEA10;Aspergillus_fumigatus_HMR_AF_270;Aspergillus_fumigatus_Z5;Aspergillus_oerlinghausenensis_CBS139183"),
            call(0.0006,"Aspergillus_fischeri_NRRL4585;Aspergillus_fumigatus_Af293;Aspergillus_fumigatus_CEA10;Aspergillus_fumigatus_HMR_AF_270;Aspergillus_fumigatus_Z5;Aspergillus_oerlinghausenensis_CBS139183"),
            call(0.0158,"Aspergillus_fumigatus_Af293;Aspergillus_fumigatus_CEA10;Aspergillus_fumigatus_HMR_AF_270;Aspergillus_fumigatus_Z5;Aspergillus_oerlinghausenensis_CBS139183"),
            call(0.0337,"Aspergillus_fumigatus_Af293;Aspergillus_fumigatus_CEA10;Aspergillus_fumigatus_HMR_AF_270;Aspergillus_fumigatus_Z5"),
            call(0.0008,"Aspergillus_fumigatus_Af293;Aspergillus_fumigatus_CEA10"),
            call(0.0007,"Aspergillus_fumigatus_HMR_AF_270;Aspergillus_fumigatus_Z5"),
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