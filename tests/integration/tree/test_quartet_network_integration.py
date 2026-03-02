import json
import sys

import pytest
from mock import patch

from phykit.phykit import Phykit


@pytest.mark.integration
class TestQuartetNetworkIntegration:
    @patch("builtins.print")
    def test_quartet_network_basic(self, mocked_print, tmp_path):
        trees = tmp_path / "trees.nwk"
        trees.write_text(
            "((A,B),(C,D));\n((A,B),(C,D));\n((A,C),(B,D));\n"
            "((A,B),(C,D));\n((A,B),(C,D));\n((A,B),(C,D));\n"
        )

        testargs = ["phykit", "quartet_network", "-t", str(trees)]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = "\n".join(str(call) for call in mocked_print.call_args_list)
        assert "Number of input trees: 6" in output
        assert "Total quartets: 1" in output

    @patch("builtins.print")
    def test_qnet_alias(self, mocked_print, tmp_path):
        trees = tmp_path / "trees.nwk"
        trees.write_text(
            "((A,B),(C,D));\n((A,B),(C,D));\n((A,B),(C,D));\n"
            "((A,B),(C,D));\n((A,B),(C,D));\n((A,B),(C,D));\n"
        )

        testargs = ["phykit", "qnet", "-t", str(trees)]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = "\n".join(str(call) for call in mocked_print.call_args_list)
        assert "Number of input trees: 6" in output

    @patch("builtins.print")
    def test_nanuq_alias(self, mocked_print, tmp_path):
        trees = tmp_path / "trees.nwk"
        trees.write_text(
            "((A,B),(C,D));\n((A,B),(C,D));\n((A,B),(C,D));\n"
            "((A,B),(C,D));\n((A,B),(C,D));\n((A,B),(C,D));\n"
        )

        testargs = ["phykit", "nanuq", "-t", str(trees)]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = "\n".join(str(call) for call in mocked_print.call_args_list)
        assert "Number of input trees: 6" in output

    @patch("builtins.print")
    def test_quartet_net_alias(self, mocked_print, tmp_path):
        trees = tmp_path / "trees.nwk"
        trees.write_text(
            "((A,B),(C,D));\n((A,B),(C,D));\n((A,B),(C,D));\n"
            "((A,B),(C,D));\n((A,B),(C,D));\n((A,B),(C,D));\n"
        )

        testargs = ["phykit", "quartet_net", "-t", str(trees)]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = "\n".join(str(call) for call in mocked_print.call_args_list)
        assert "Number of input trees: 6" in output

    @patch("builtins.print")
    def test_json_output(self, mocked_print, tmp_path):
        trees = tmp_path / "trees.nwk"
        trees.write_text(
            "((A,B),(C,D));\n((A,B),(C,D));\n((A,C),(B,D));\n"
            "((A,B),(C,D));\n((A,B),(C,D));\n((A,B),(C,D));\n"
        )

        testargs = ["phykit", "quartet_network", "-t", str(trees), "--json"]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["input_tree_count"] == 6
        assert payload["total_quartets"] == 1
        assert len(payload["quartets"]) == 1
        assert "p_star" in payload["quartets"][0]
        assert "p_tree" in payload["quartets"][0]
        assert payload["beta"] == 0.95

    @patch("builtins.print")
    def test_custom_alpha(self, mocked_print, tmp_path):
        trees = tmp_path / "trees.nwk"
        trees.write_text(
            "((A,B),(C,D));\n((A,B),(C,D));\n((A,C),(B,D));\n"
            "((A,B),(C,D));\n((A,B),(C,D));\n((A,B),(C,D));\n"
        )

        testargs = [
            "phykit", "quartet_network", "-t", str(trees),
            "--alpha", "0.01", "--json"
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["alpha"] == 0.01

    @patch("builtins.print")
    def test_custom_beta(self, mocked_print, tmp_path):
        trees = tmp_path / "trees.nwk"
        trees.write_text(
            "((A,B),(C,D));\n((A,B),(C,D));\n((A,C),(B,D));\n"
            "((A,B),(C,D));\n((A,B),(C,D));\n((A,B),(C,D));\n"
        )

        testargs = [
            "phykit", "quartet_network", "-t", str(trees),
            "--beta", "0.90", "--json"
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["beta"] == 0.90

    @patch("builtins.print")
    def test_plot_output_creates_file(self, mocked_print, tmp_path):
        trees = tmp_path / "trees.nwk"
        trees.write_text(
            "((A,B),(C,D));\n((A,B),(C,D));\n((A,C),(B,D));\n"
            "((A,B),(C,D));\n((A,B),(C,D));\n((A,B),(C,D));\n"
        )
        plot_file = tmp_path / "qnet.png"

        testargs = [
            "phykit", "quartet_network", "-t", str(trees),
            "--plot-output", str(plot_file)
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert plot_file.exists()

    @patch("builtins.print")
    def test_sample_file(self, mocked_print):
        testargs = [
            "phykit", "qnet", "-t",
            "tests/sample_files/gene_trees_for_network.nwk",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = "\n".join(str(call) for call in mocked_print.call_args_list)
        assert "Number of input trees:" in output
        assert "Total quartets:" in output
