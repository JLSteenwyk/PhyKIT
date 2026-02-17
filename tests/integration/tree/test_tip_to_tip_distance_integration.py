from mock import patch, call
from pathlib import Path
import pytest
import sys
import json

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
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_tip_to_tip_distance_bad_file_path(self, mocked_print):
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_tip_to_tip_distance_json(self, mocked_print):
        testargs = [
            "phykit",
            "tip_to_tip_distance",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tre_rooted.tree",
            "Aspergillus_fumigatus_HMR_AF_270",
            "Aspergillus_fischeri_NRRL181.GCF_000149645.1_ASM14964v1",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload == {
            "taxon_a": "Aspergillus_fumigatus_HMR_AF_270",
            "taxon_b": "Aspergillus_fischeri_NRRL181.GCF_000149645.1_ASM14964v1",
            "tip_to_tip_distance": 0.0539,
        }

    @patch("builtins.print")
    def test_tip_to_tip_distance_all_pairs(self, mocked_print):
        testargs = [
            "phykit",
            "tip_to_tip_distance",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tre_rooted.tree",
            "--all-pairs",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert len(mocked_print.mock_calls) == 45
        assert any(
            call_args.args[0] == (
                "Aspergillus_fumigatus_Af293\tAspergillus_fumigatus_CEA10\t0.0021"
            )
            for call_args in mocked_print.mock_calls
        )

    @patch("builtins.print")
    def test_tip_to_tip_distance_all_pairs_json(self, mocked_print):
        testargs = [
            "phykit",
            "tip_to_tip_distance",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tre_rooted.tree",
            "--all-pairs",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["all_pairs"] is True
        assert payload["rows"][0] == payload["pairs"][0]
        assert len(payload["rows"]) == 45

    @patch("phykit.services.tree.tip_to_tip_distance.TipToTipDistance._plot_tip_distance_heatmap")
    @patch("builtins.print")
    def test_tip_to_tip_distance_all_pairs_plot(self, mocked_print, mocked_plot):
        testargs = [
            "phykit",
            "tip_to_tip_distance",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tre_rooted.tree",
            "--all-pairs",
            "--plot",
            "--plot-output",
            "t2t_heatmap_test.png",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_plot.called
        assert any(
            call_args.args[0] == "Saved tip-to-tip heatmap: t2t_heatmap_test.png"
            for call_args in mocked_print.mock_calls
        )

    @patch("phykit.services.tree.tip_to_tip_distance.TipToTipDistance._plot_tip_distance_heatmap")
    @patch("builtins.print")
    def test_tip_to_tip_distance_all_pairs_plot_json(self, mocked_print, mocked_plot):
        testargs = [
            "phykit",
            "tip_to_tip_distance",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tre_rooted.tree",
            "--all-pairs",
            "--plot",
            "--plot-output",
            "t2t_heatmap_test_json.png",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_plot.called
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["plot_output"] == "t2t_heatmap_test_json.png"

    @patch("builtins.print")
    def test_tip_to_tip_distance_plot_requires_all_pairs(self, mocked_print):
        testargs = [
            "phykit",
            "tip_to_tip_distance",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tre_rooted.tree",
            "Aspergillus_fumigatus_Af293",
            "Aspergillus_fischeri_IBT_3007",
            "--plot",
        ]
        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 2
