import pytest
import sys
from mock import patch, call
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestOccupancyPerTaxon(object):
    @patch("builtins.print")
    def test_occupancy_per_taxon0(self, mocked_print):
        expected_result_0 = "1\t0.8333"
        expected_result_1 = "2\t0.6667"
        expected_result_2 = "3\t0.6667"
        expected_result_3 = "4\t0.8333"
        expected_result_4 = "5\t0.6667"
        testargs = [
            "phykit",
            "occupancy_per_taxon",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
            call(expected_result_2),
            call(expected_result_3),
            call(expected_result_4),
        ]

    @patch("builtins.print")
    def test_occupancy_per_taxon_alias(self, mocked_print):
        expected_result_0 = "1\t0.8333"
        expected_result_1 = "2\t0.6667"
        expected_result_2 = "3\t0.6667"
        expected_result_3 = "4\t0.8333"
        expected_result_4 = "5\t0.6667"
        testargs = [
            "phykit",
            "occ_tax",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
            call(expected_result_2),
            call(expected_result_3),
            call(expected_result_4),
        ]

    @patch("builtins.print")
    def test_occupancy_per_taxon_incorrect_input_file(self, mocked_print):
        testargs = [
            "phykit",
            "occupancy_per_taxon",
            f"{here.parent.parent.parent}/sample_files/test_alignment_1.f",
        ]
        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2
