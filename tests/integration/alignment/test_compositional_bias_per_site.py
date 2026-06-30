import pytest
import sys
import json
from mock import patch, call
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)


EXPECTED_SIMPLE_CBPS = (
    "1\t0.0\tnan\tnan\n"
    "2\t0.0\t1.0\t1.0\n"
    "3\t0.2\t1.0\t0.6547\n"
    "4\t0.0\tnan\tnan\n"
    "5\t0.2\t1.0\t0.6547\n"
    "6\t0.0\t1.0\t1.0"
)


@pytest.mark.integration
class TestCompositionalBiasPerSite(object):
    @patch("builtins.print")
    def test_compositional_bias_per_site(self, mocked_print):
        testargs = [
            "phykit",
            "compositional_bias_per_site",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(EXPECTED_SIMPLE_CBPS)]

    @patch("builtins.print")
    def test_compositional_bias_per_site_alias0(self, mocked_print):
        testargs = [
            "phykit",
            "comp_bias_per_site",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(EXPECTED_SIMPLE_CBPS)]

    @patch("builtins.print")
    def test_compositional_bias_per_site_alias1(self, mocked_print):
        testargs = [
            "phykit",
            "cbps",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(EXPECTED_SIMPLE_CBPS)]

    @patch("builtins.print")
    def test_compositional_bias_per_site_json(self, mocked_print):
        testargs = [
            "phykit",
            "compositional_bias_per_site",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["rows"][0] == payload["sites"][0]
        assert payload["sites"][0] == {
            "chi_square": 0.0,
            "p_value": None,
            "p_value_corrected": None,
            "site": 1,
        }
        assert payload["sites"][2] == {
            "chi_square": 0.2,
            "p_value": 0.6547,
            "p_value_corrected": 1.0,
            "site": 3,
        }

    @patch("phykit.services.alignment.compositional_bias_per_site.CompositionalBiasPerSite._plot_compositional_bias_manhattan")
    @patch("builtins.print")
    def test_compositional_bias_per_site_plot(self, mocked_print, mocked_plot):
        testargs = [
            "phykit",
            "compositional_bias_per_site",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "--plot",
            "--plot-output",
            "cbps_test_plot.png",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_plot.called
        assert any(
            call.args[0] == "Saved compositional bias plot: cbps_test_plot.png"
            for call in mocked_print.mock_calls
        )

    @patch("phykit.services.alignment.compositional_bias_per_site.CompositionalBiasPerSite._plot_compositional_bias_manhattan")
    @patch("builtins.print")
    def test_compositional_bias_per_site_plot_json(self, mocked_print, mocked_plot):
        testargs = [
            "phykit",
            "compositional_bias_per_site",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "--plot",
            "--plot-output",
            "cbps_test_plot_json.png",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_plot.called
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["plot_output"] == "cbps_test_plot_json.png"
