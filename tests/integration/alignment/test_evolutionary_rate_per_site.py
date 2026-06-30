import pytest
import sys
import json
from mock import patch, call
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)


EXPECTED_SIMPLE_ERPS = (
    "1\t0.0\n"
    "2\t0.5\n"
    "3\t0.48\n"
    "4\t0.0\n"
    "5\t0.48\n"
    "6\t0.5"
)


@pytest.mark.integration
class TestEvolutionaryRatePerSite(object):
    @patch("builtins.print")
    def test_evolutionary_rate_per_site(self, mocked_print):
        testargs = [
            "phykit",
            "evolutionary_rate_per_site",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(EXPECTED_SIMPLE_ERPS)]

    @patch("builtins.print")
    def test_evolutionary_rate_per_site_alias0(self, mocked_print):
        testargs = [
            "phykit",
            "evo_rate_per_site",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(EXPECTED_SIMPLE_ERPS)]

    @patch("builtins.print")
    def test_evolutionary_rate_per_site_alias1(self, mocked_print):
        testargs = [
            "phykit",
            "erps",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(EXPECTED_SIMPLE_ERPS)]

    @patch("builtins.print")
    def test_evolutionary_rate_per_site_json(self, mocked_print):
        testargs = [
            "phykit",
            "evolutionary_rate_per_site",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["rows"][0] == payload["sites"][0]
        assert payload["sites"][0] == {"evolutionary_rate": 0.0, "site": 1}
        assert payload["sites"][2] == {"evolutionary_rate": 0.48, "site": 3}

    @patch("phykit.services.alignment.evolutionary_rate_per_site.EvolutionaryRatePerSite._plot_evolutionary_rate_per_site")
    @patch("builtins.print")
    def test_evolutionary_rate_per_site_plot(self, mocked_print, mocked_plot):
        testargs = [
            "phykit",
            "evolutionary_rate_per_site",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "--plot",
            "--plot-output",
            "erps_test_plot.png",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_plot.called
        assert any(
            call.args[0] == "Saved evolutionary-rate plot: erps_test_plot.png"
            for call in mocked_print.mock_calls
        )

    @patch("phykit.services.alignment.evolutionary_rate_per_site.EvolutionaryRatePerSite._plot_evolutionary_rate_per_site")
    @patch("builtins.print")
    def test_evolutionary_rate_per_site_plot_json(self, mocked_print, mocked_plot):
        testargs = [
            "phykit",
            "evolutionary_rate_per_site",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "--plot",
            "--plot-output",
            "erps_test_plot_json.png",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_plot.called
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["plot_output"] == "erps_test_plot_json.png"
