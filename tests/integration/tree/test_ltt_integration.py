import json
import sys

import pytest
from mock import patch

from phykit.phykit import Phykit


BALANCED_TREE = "tests/sample_files/ltt_test_balanced.tre"
LADDER_TREE = "tests/sample_files/ltt_test_ladder.tre"
RECENT_BURST = "tests/sample_files/ltt_test_recent_burst.tre"


def _output(mocked_print):
    """Collect all print() calls into one string."""
    return "\n".join(str(call) for call in mocked_print.call_args_list)


@pytest.mark.integration
class TestLTTIntegration:
    @patch("builtins.print")
    def test_basic_ltt(self, mocked_print):
        testargs = ["phykit", "ltt", "-t", BALANCED_TREE]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = _output(mocked_print)
        assert "-1.4142" in output
        assert "0.1573" in output

    @patch("builtins.print")
    def test_gamma_stat_alias(self, mocked_print):
        testargs = ["phykit", "gamma_stat", "-t", BALANCED_TREE]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = _output(mocked_print)
        assert "-1.4142" in output

    @patch("builtins.print")
    def test_gamma_alias(self, mocked_print):
        testargs = ["phykit", "gamma", "-t", LADDER_TREE]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = _output(mocked_print)
        assert "-0.7143" in output

    @patch("builtins.print")
    def test_verbose_output(self, mocked_print):
        testargs = ["phykit", "ltt", "-t", BALANCED_TREE, "-v"]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = _output(mocked_print)
        assert "-1.4142" in output
        assert "Branching times" in output
        assert "Lineage-through-time" in output

    @patch("builtins.print")
    def test_json_output(self, mocked_print):
        testargs = ["phykit", "ltt", "-t", RECENT_BURST, "--json"]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert "gamma" in payload
        assert "p_value" in payload
        assert "ltt" in payload
        assert "branching_times" in payload
        assert abs(payload["gamma"] - 2.2825) < 0.001

    @patch("builtins.print")
    def test_plot_output(self, mocked_print, tmp_path):
        plot_file = tmp_path / "ltt_plot.png"
        testargs = [
            "phykit", "ltt", "-t", BALANCED_TREE,
            "--plot-output", str(plot_file),
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert plot_file.exists()
        assert plot_file.stat().st_size > 0

    @patch("builtins.print")
    def test_recent_burst_positive_gamma(self, mocked_print):
        testargs = ["phykit", "ltt", "-t", RECENT_BURST]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = _output(mocked_print)
        assert "2.2825" in output

    def test_too_few_tips(self, tmp_path):
        tree_file = tmp_path / "tiny.tre"
        tree_file.write_text("(A:1,B:1);")
        testargs = ["phykit", "ltt", "-t", str(tree_file)]
        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit):
                Phykit()
