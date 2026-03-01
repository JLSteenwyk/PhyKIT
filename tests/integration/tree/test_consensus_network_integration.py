import json
import sys

import pytest
from mock import patch

from phykit.phykit import Phykit


@pytest.mark.integration
class TestConsensusNetworkIntegration:
    @patch("builtins.print")
    def test_consensus_network_basic(self, mocked_print, tmp_path):
        trees = tmp_path / "trees.nwk"
        trees.write_text("((A,B),(C,D));\n((A,B),(C,D));\n((A,C),(B,D));\n")

        testargs = ["phykit", "consensus_network", "-t", str(trees)]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = mocked_print.call_args_list
        full_output = "\n".join(str(call) for call in output)
        assert "Number of input trees: 3" in full_output

    @patch("builtins.print")
    def test_consnet_alias(self, mocked_print, tmp_path):
        trees = tmp_path / "trees.nwk"
        trees.write_text("((A,B),(C,D));\n((A,B),(C,D));\n")

        testargs = ["phykit", "consnet", "-t", str(trees)]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = "\n".join(str(call) for call in mocked_print.call_args_list)
        assert "Number of input trees: 2" in output

    @patch("builtins.print")
    def test_splitnet_alias(self, mocked_print, tmp_path):
        trees = tmp_path / "trees.nwk"
        trees.write_text("((A,B),(C,D));\n((A,B),(C,D));\n")

        testargs = ["phykit", "splitnet", "-t", str(trees)]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = "\n".join(str(call) for call in mocked_print.call_args_list)
        assert "Number of input trees: 2" in output

    @patch("builtins.print")
    def test_splits_network_alias(self, mocked_print, tmp_path):
        trees = tmp_path / "trees.nwk"
        trees.write_text("((A,B),(C,D));\n((A,B),(C,D));\n")

        testargs = ["phykit", "splits_network", "-t", str(trees)]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = "\n".join(str(call) for call in mocked_print.call_args_list)
        assert "Number of input trees: 2" in output

    @patch("builtins.print")
    def test_json_output(self, mocked_print, tmp_path):
        trees = tmp_path / "trees.nwk"
        trees.write_text("((A,B),(C,D));\n((A,B),(C,D));\n")

        testargs = ["phykit", "consensus_network", "-t", str(trees), "--json"]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["input_tree_count"] == 2
        assert payload["taxa_count"] == 4
        assert len(payload["splits"]) > 0

    @patch("builtins.print")
    def test_custom_threshold(self, mocked_print, tmp_path):
        trees = tmp_path / "trees.nwk"
        trees.write_text(
            "((A,B),(C,D));\n((A,B),(C,D));\n((A,C),(B,D));\n"
        )

        testargs = [
            "phykit", "consensus_network", "-t", str(trees),
            "--threshold", "0.7", "--json"
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        # Only splits with freq >= 0.7 should be included
        for s in payload["splits"]:
            assert s["frequency"] >= 0.7

    @patch("builtins.print")
    def test_plot_output_creates_file(self, mocked_print, tmp_path):
        trees = tmp_path / "trees.nwk"
        trees.write_text("((A,B),(C,D));\n((A,B),(C,D));\n")
        plot_file = tmp_path / "network.png"

        testargs = [
            "phykit", "consensus_network", "-t", str(trees),
            "--plot-output", str(plot_file)
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert plot_file.exists()

    @patch("builtins.print")
    def test_missing_taxa_shared(self, mocked_print, tmp_path):
        trees = tmp_path / "trees.nwk"
        trees.write_text("((A,B),(C,D));\n((A,B),(C,E));\n")

        testargs = [
            "phykit", "consnet", "-t", str(trees),
            "--missing-taxa", "shared", "--json"
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["pruned_to_shared_taxa"] is True
        assert payload["taxa_count"] == 3

    @patch("builtins.print")
    def test_sample_file(self, mocked_print):
        testargs = [
            "phykit", "consnet", "-t",
            "tests/sample_files/gene_trees_for_network.nwk",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = "\n".join(str(call) for call in mocked_print.call_args_list)
        assert "Number of input trees:" in output
