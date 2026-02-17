import pytest
import sys
import json
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestRCVT(object):
    @patch("builtins.print")
    def test_rcvt0(self, mocked_print):
        expected_result_0 = dedent(
            """1\t0.056"""
        )
        expected_result_1 = dedent(
            """2\t0.04"""
        )
        expected_result_2 = dedent(
            """3\t0.04"""  
        )
        expected_result_3 = dedent(
            """4\t0.056"""
        )
        expected_result_4 = dedent(
            """5\t0.1"""  
        )
        testargs = [
            "phykit",
            "rcvt",
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
    def test_rcvt1(self, mocked_print):
        expected_result_0 = dedent(
            """1\t0.0781"""
        )
        expected_result_1 = dedent(
            """2\t0.0556"""
        )
        expected_result_2 = dedent(
            """3\t0.0278"""  
        )
        expected_result_3 = dedent(
            """4\t0.1094"""
        )
        testargs = [
            "phykit",
            "rcvt",
            f"{here.parent.parent.parent}/sample_files/test_alignment_0.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
            call(expected_result_2),
            call(expected_result_3),
        ]

    @patch("builtins.print")
    def test_rcvt2(self, mocked_print):
        expected_result_0 = dedent(
            """1\t0.0417"""
        )
        expected_result_1 = dedent(
            """2\t0.0417"""
        )
        expected_result_2 = dedent(
            """3\t0.0417"""  
        )
        expected_result_3 = dedent(
            """4\t0.125"""
        )
        testargs = [
            "phykit",
            "rcvt",
            f"{here.parent.parent.parent}/sample_files/test_alignment_1.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
            call(expected_result_2),
            call(expected_result_3),
        ]

    @patch("builtins.print")
    def test_rcvt_alias(self, mocked_print):
        expected_result_0 = dedent(
            """1\t0.0417"""
        )
        expected_result_1 = dedent(
            """2\t0.0417"""
        )
        expected_result_2 = dedent(
            """3\t0.0417"""  
        )
        expected_result_3 = dedent(
            """4\t0.125"""
        )
        testargs = [
            "phykit",
            "relative_composition_variability_taxon",
            f"{here.parent.parent.parent}/sample_files/test_alignment_1.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
            call(expected_result_2),
            call(expected_result_3),
        ]

    @patch("builtins.print")
    def test_rcvt_incorrect_input_file(self, mocked_print):
        testargs = [
            "phykit",
            "rcvt",
            f"{here.parent.parent.parent}/sample_files/does_not_exist",
        ]
        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_rcvt_json(self, mocked_print):
        testargs = [
            "phykit",
            "rcvt",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["rows"][0] == payload["taxa"][0]
        assert payload["taxa"][0] == {"taxon": "1", "rcvt": 0.056}

    @patch("phykit.services.alignment.rcvt.RelativeCompositionVariabilityTaxon._plot_rcvt")
    @patch("builtins.print")
    def test_rcvt_plot(self, mocked_print, mocked_plot):
        testargs = [
            "phykit",
            "rcvt",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "--plot",
            "--plot-output",
            "rcvt_test_plot.png",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_plot.called
        assert any(
            call.args[0] == "Saved RCVT plot: rcvt_test_plot.png"
            for call in mocked_print.mock_calls
        )

    @patch("phykit.services.alignment.rcvt.RelativeCompositionVariabilityTaxon._plot_rcvt")
    @patch("builtins.print")
    def test_rcvt_plot_json(self, mocked_print, mocked_plot):
        testargs = [
            "phykit",
            "rcvt",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "--plot",
            "--plot-output",
            "rcvt_test_plot_json.png",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_plot.called
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["plot_output"] == "rcvt_test_plot_json.png"
