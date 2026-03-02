import json
import sys

import pytest
from mock import patch

from phykit.phykit import Phykit


TREE_SIMPLE = "tests/sample_files/tree_simple.tre"
TRAITS_FILE = "tests/sample_files/tree_simple_threshold_traits.tsv"


def _output(mocked_print):
    """Collect all print() calls into one string."""
    return "\n".join(str(call) for call in mocked_print.call_args_list)


@pytest.mark.integration
class TestThresholdModelIntegration:
    @patch("builtins.print")
    def test_threshold_model_basic(self, mocked_print):
        testargs = [
            "phykit", "threshold_model",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
            "--traits", "habitat,body_mass",
            "--types", "discrete,continuous",
            "--ngen", "2000",
            "--sample", "10",
            "--seed", "42",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = _output(mocked_print)
        assert "Posterior correlation (r):" in output
        assert "Trait 1: habitat" in output

    @patch("builtins.print")
    def test_threshold_alias(self, mocked_print):
        testargs = [
            "phykit", "threshold",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
            "--traits", "habitat,body_mass",
            "--types", "discrete,continuous",
            "--ngen", "1000",
            "--sample", "10",
            "--seed", "42",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = _output(mocked_print)
        assert "Posterior correlation (r):" in output

    @patch("builtins.print")
    def test_thresh_alias(self, mocked_print):
        testargs = [
            "phykit", "thresh",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
            "--traits", "habitat,body_mass",
            "--types", "discrete,continuous",
            "--ngen", "1000",
            "--sample", "10",
            "--seed", "42",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = _output(mocked_print)
        assert "Posterior correlation (r):" in output

    @patch("builtins.print")
    def test_threshbayes_alias(self, mocked_print):
        testargs = [
            "phykit", "threshbayes",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
            "--traits", "habitat,body_mass",
            "--types", "discrete,continuous",
            "--ngen", "1000",
            "--sample", "10",
            "--seed", "42",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = _output(mocked_print)
        assert "Posterior correlation (r):" in output

    @patch("builtins.print")
    def test_json_output(self, mocked_print):
        testargs = [
            "phykit", "threshold_model",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
            "--traits", "habitat,body_mass",
            "--types", "discrete,continuous",
            "--ngen", "2000",
            "--sample", "10",
            "--seed", "42",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert "metadata" in payload
        assert "summary" in payload
        assert "posterior_samples" in payload
        assert payload["metadata"]["trait1"] == "habitat"
        assert "r" in payload["summary"]

    @patch("builtins.print")
    def test_plot_output(self, mocked_print, tmp_path):
        plot_file = tmp_path / "trace.png"
        testargs = [
            "phykit", "threshold_model",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
            "--traits", "habitat,body_mass",
            "--types", "discrete,continuous",
            "--ngen", "2000",
            "--sample", "10",
            "--seed", "42",
            "--plot", str(plot_file),
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert plot_file.exists()
        assert plot_file.stat().st_size > 0
