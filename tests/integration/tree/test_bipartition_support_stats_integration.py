import pytest
import sys
import json
from mock import patch
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)

SUMMARY_OUTPUT = (
    "mean: 95.7143\n"
    "median: 100.0\n"
    "25th percentile: 92.5\n"
    "75th percentile: 100.0\n"
    "minimum: 85.0\n"
    "maximum: 100.0\n"
    "standard deviation: 7.3193\n"
    "variance: 53.5714"
)


@pytest.mark.integration
class TestBipartitionSupportStats(object):
    @patch("builtins.print")
    def test_bipartition_support_stats(self, mocked_print):
        testargs = [
            "phykit",
            "bipartition_support_stats",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        mocked_print.assert_called_once_with(SUMMARY_OUTPUT)

    @patch("builtins.print")
    def test_bipartition_support_stats_verbose(self, mocked_print):
        testargs = [
            "phykit",
            "bipartition_support_stats",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            "-v"
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.call_args.args[0].splitlines() == [
            "85 Aspergillus_fischeri_IBT_3007;Aspergillus_fischeri_NRRL181.GCF_000149645.1_ASM14964v1;Aspergillus_fischeri_NRRL4585;Aspergillus_fumigatus_Af293;Aspergillus_fumigatus_CEA10;Aspergillus_fumigatus_HMR_AF_270;Aspergillus_fumigatus_Z5;Aspergillus_oerlinghausenensis_CBS139183",
            "85 Aspergillus_fischeri_NRRL181.GCF_000149645.1_ASM14964v1;Aspergillus_fischeri_NRRL4585;Aspergillus_fumigatus_Af293;Aspergillus_fumigatus_CEA10;Aspergillus_fumigatus_HMR_AF_270;Aspergillus_fumigatus_Z5;Aspergillus_oerlinghausenensis_CBS139183",
            "100 Aspergillus_fischeri_NRRL4585;Aspergillus_fumigatus_Af293;Aspergillus_fumigatus_CEA10;Aspergillus_fumigatus_HMR_AF_270;Aspergillus_fumigatus_Z5;Aspergillus_oerlinghausenensis_CBS139183",
            "100 Aspergillus_fumigatus_Af293;Aspergillus_fumigatus_CEA10;Aspergillus_fumigatus_HMR_AF_270;Aspergillus_fumigatus_Z5;Aspergillus_oerlinghausenensis_CBS139183",
            "100 Aspergillus_fumigatus_Af293;Aspergillus_fumigatus_CEA10;Aspergillus_fumigatus_HMR_AF_270;Aspergillus_fumigatus_Z5",
            "100 Aspergillus_fumigatus_Af293;Aspergillus_fumigatus_CEA10",
            "100 Aspergillus_fumigatus_HMR_AF_270;Aspergillus_fumigatus_Z5",
        ]

    @patch("builtins.print")
    def test_bipartition_support_incorrect_input_file(self, mocked_print):
        testargs = [
            "phykit",
            "bipartition_support_stats",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tr"
        ]

        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_bipartition_support_stats_alias(self, mocked_print):
        testargs = [
            "phykit",
            "bss",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        mocked_print.assert_called_once_with(SUMMARY_OUTPUT)

    @patch("builtins.print")
    def test_bipartition_support_stats_thresholds(self, mocked_print):
        testargs = [
            "phykit",
            "bipartition_support_stats",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            "--thresholds",
            "90,100",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.call_args_list[-2].args[0] == SUMMARY_OUTPUT
        assert mocked_print.call_args.args[0] == "\n".join(
            [
                "below 90.0: 2 (28.5714%)",
                "below 100.0: 2 (28.5714%)",
            ]
        )

    @patch("builtins.print")
    def test_bipartition_support_stats_json(self, mocked_print):
        testargs = [
            "phykit",
            "bipartition_support_stats",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            "--thresholds",
            "90",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["verbose"] is False
        assert round(payload["summary"]["mean"], 4) == 95.7143
        assert payload["thresholds"] == [
            {"count_below": 2, "fraction_below": 2 / 7, "threshold": 90.0}
        ]
