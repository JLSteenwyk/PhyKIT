import pytest
import sys
import json
from mock import patch, call
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestCompositionPerTaxon(object):
    @patch("builtins.print")
    def test_composition_per_taxon0(self, mocked_print):
        expected_result_0 = "1\tA:0.4;C:0.0;G:0.2;T:0.4"
        expected_result_1 = "2\tA:0.5;C:0.0;G:0.25;T:0.25"
        expected_result_2 = "3\tA:0.5;C:0.0;G:0.25;T:0.25"
        expected_result_3 = "4\tA:0.6;C:0.0;G:0.2;T:0.2"
        expected_result_4 = "5\tA:0.5;C:0.25;G:0.0;T:0.25"
        testargs = [
            "phykit",
            "composition_per_taxon",
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
    def test_composition_per_taxon_alias(self, mocked_print):
        expected_result_0 = "1\tA:0.4;C:0.0;G:0.2;T:0.4"
        expected_result_1 = "2\tA:0.5;C:0.0;G:0.25;T:0.25"
        expected_result_2 = "3\tA:0.5;C:0.0;G:0.25;T:0.25"
        expected_result_3 = "4\tA:0.6;C:0.0;G:0.2;T:0.2"
        expected_result_4 = "5\tA:0.5;C:0.25;G:0.0;T:0.25"
        testargs = [
            "phykit",
            "comp_tax",
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
    def test_composition_per_taxon_incorrect_input_file(self, mocked_print):
        testargs = [
            "phykit",
            "composition_per_taxon",
            f"{here.parent.parent.parent}/sample_files/test_alignment_1.f",
        ]
        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_composition_per_taxon_json(self, mocked_print):
        testargs = [
            "phykit",
            "composition_per_taxon",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["symbols"] == ["A", "C", "G", "T"]
        assert payload["rows"][0] == payload["taxa"][0]
        assert payload["taxa"][0] == {
            "taxon": "1",
            "composition": {"A": 0.4, "C": 0.0, "G": 0.2, "T": 0.4},
        }
