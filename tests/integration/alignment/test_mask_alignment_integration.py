import pytest
import sys
import json
from mock import patch, call
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestMaskAlignment(object):
    @patch("builtins.print")
    def test_mask_alignment_with_gap_and_occupancy_thresholds(self, mocked_print):
        testargs = [
            "phykit",
            "mask_alignment",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-g",
            "0.3",
            "-o",
            "0.8",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(">1\nAGAT"),
            call(">2\nAGAT"),
            call(">3\nAGTA"),
            call(">4\nAATA"),
            call(">5\nAAT-"),
        ]

    @patch("builtins.print")
    def test_mask_alignment_alias(self, mocked_print):
        testargs = [
            "phykit",
            "mask",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-e",
            "0.5",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(">1\nAT"),
            call(">2\nA-"),
            call(">3\nA-"),
            call(">4\nA-"),
            call(">5\nA-"),
        ]

    @patch("builtins.print")
    def test_mask_alignment_incorrect_input_file(self, mocked_print):
        testargs = [
            "phykit",
            "mask_alignment",
            f"{here.parent.parent.parent}/sample_files/test_alignment_1.f",
        ]
        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_mask_alignment_json(self, mocked_print):
        testargs = [
            "phykit",
            "mask_alignment",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-g",
            "0.3",
            "-o",
            "0.8",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["kept_sites"] == 4
        assert payload["total_sites"] == 6
        assert payload["thresholds"]["max_gap"] == 0.3
        assert payload["rows"][0] == payload["taxa"][0]
        assert payload["taxa"][0] == {"taxon": "1", "sequence": "AGAT"}
