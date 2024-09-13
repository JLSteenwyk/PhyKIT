import pytest
import sys
from mock import patch, call
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)


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

        assert mocked_print.mock_calls == [
            call("mean: 95.7143"),
            call("median: 100"),
            call("25th percentile: 92.5"),
            call("75th percentile: 100.0"),
            call("minimum: 85"),
            call("maximum: 100"),
            call("standard deviation: 7.3193"),
            call("variance: 53.5714")
        ]

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

        assert mocked_print.mock_calls == [
            call(85, "Aspergillus_fischeri_IBT_3007;Aspergillus_fischeri_NRRL181.GCF_000149645.1_ASM14964v1;Aspergillus_fischeri_NRRL4585;Aspergillus_fumigatus_Af293;Aspergillus_fumigatus_CEA10;Aspergillus_fumigatus_HMR_AF_270;Aspergillus_fumigatus_Z5;Aspergillus_oerlinghausenensis_CBS139183"),
            call(85, "Aspergillus_fischeri_NRRL181.GCF_000149645.1_ASM14964v1;Aspergillus_fischeri_NRRL4585;Aspergillus_fumigatus_Af293;Aspergillus_fumigatus_CEA10;Aspergillus_fumigatus_HMR_AF_270;Aspergillus_fumigatus_Z5;Aspergillus_oerlinghausenensis_CBS139183"),
            call(100, "Aspergillus_fischeri_NRRL4585;Aspergillus_fumigatus_Af293;Aspergillus_fumigatus_CEA10;Aspergillus_fumigatus_HMR_AF_270;Aspergillus_fumigatus_Z5;Aspergillus_oerlinghausenensis_CBS139183"),
            call(100, "Aspergillus_fumigatus_Af293;Aspergillus_fumigatus_CEA10;Aspergillus_fumigatus_HMR_AF_270;Aspergillus_fumigatus_Z5;Aspergillus_oerlinghausenensis_CBS139183"),
            call(100, "Aspergillus_fumigatus_Af293;Aspergillus_fumigatus_CEA10;Aspergillus_fumigatus_HMR_AF_270;Aspergillus_fumigatus_Z5"),
            call(100, "Aspergillus_fumigatus_Af293;Aspergillus_fumigatus_CEA10"),
            call(100, "Aspergillus_fumigatus_HMR_AF_270;Aspergillus_fumigatus_Z5"),
        ]

    @patch("builtins.print")
    def test_bipartition_support_incorrect_input_file(self, mocked_print):
        testargs = [
            "phykit",
            "bipartition_support_stats",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tr"
        ]

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
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

        assert mocked_print.mock_calls == [
            call("mean: 95.7143"),
            call("median: 100"),
            call("25th percentile: 92.5"),
            call("75th percentile: 100.0"),
            call("minimum: 85"),
            call("maximum: 100"),
            call("standard deviation: 7.3193"),
            call("variance: 53.5714")
        ]
