from mock import patch, call
from pathlib import Path
import pytest
import sys

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestTipToTipDistance(object):
    @patch("builtins.print")
    def test_tip_to_tip_distance_0(self, mocked_print):
        expected_result = 0.0532
        testargs = [
            "phykit",
            "tip_to_tip_distance",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tre_rooted.tree",
            "Aspergillus_fumigatus_Af293",
            "Aspergillus_fischeri_IBT_3007",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_tip_to_tip_distance_1(self, mocked_print):
        expected_result = 0.0539
        testargs = [
            "phykit",
            "tip_to_tip_distance",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tre_rooted.tree",
            "Aspergillus_fumigatus_HMR_AF_270",
            "Aspergillus_fischeri_NRRL181.GCF_000149645.1_ASM14964v1",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_tip_to_tip_distance_alias_0(self, mocked_print):
        expected_result = 0.0539
        testargs = [
            "phykit",
            "t2t_dist",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tre_rooted.tree",
            "Aspergillus_fumigatus_HMR_AF_270",
            "Aspergillus_fischeri_NRRL181.GCF_000149645.1_ASM14964v1",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_tip_to_tip_distance_alias_1(self, mocked_print):
        expected_result = 0.0539
        testargs = [
            "phykit",
            "t2t",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tre_rooted.tree",
            "Aspergillus_fumigatus_HMR_AF_270",
            "Aspergillus_fischeri_NRRL181.GCF_000149645.1_ASM14964v1",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_tip_to_tip_distance_tip_not_in_tree(self, mocked_print):
        testargs = [
            "phykit",
            "t2t",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tre_rooted.tree",
            "Aspergillus_fumigatus_HMR_AF_270",
            "not_in_tree",
        ]
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_tip_to_tip_distance_bad_file_path(self, mocked_print):
        testargs = [
            "phykit",
            "t2t",
            f"{here.parent.parent.parent}/sample_files/no_tree",
            "Aspergillus_fumigatus_HMR_AF_270",
            "Aspergillus_fischeri_NRRL181.GCF_000149645.1_ASM14964v1",
        ]
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2
