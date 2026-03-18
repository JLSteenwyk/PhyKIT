import json
import sys

import pytest
from mock import patch

from phykit.phykit import Phykit


@pytest.mark.integration
class TestTreeSpaceIntegration:
    @patch("builtins.print")
    def test_basic_invocation(self, mocked_print, tmp_path):
        output = tmp_path / "tree_space.png"
        testargs = [
            "phykit", "tree_space",
            "-t", "tests/sample_files/gene_trees_simple.nwk",
            "-o", str(output),
            "--seed", "42",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output_text = "\n".join(str(call) for call in mocked_print.call_args_list)
        assert "Tree Space Analysis" in output_text
        assert "Gene trees: 10" in output_text
        assert output.exists()

    @patch("builtins.print")
    def test_alias_tspace(self, mocked_print, tmp_path):
        output = tmp_path / "tree_space.png"
        testargs = [
            "phykit", "tspace",
            "-t", "tests/sample_files/gene_trees_simple.nwk",
            "-o", str(output),
            "--seed", "42",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output_text = "\n".join(str(call) for call in mocked_print.call_args_list)
        assert "Tree Space Analysis" in output_text

    @patch("builtins.print")
    def test_alias_tree_landscape(self, mocked_print, tmp_path):
        output = tmp_path / "tree_space.png"
        testargs = [
            "phykit", "tree_landscape",
            "-t", "tests/sample_files/gene_trees_simple.nwk",
            "-o", str(output),
            "--seed", "42",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output_text = "\n".join(str(call) for call in mocked_print.call_args_list)
        assert "Tree Space Analysis" in output_text

    @patch("builtins.print")
    def test_kf_metric(self, mocked_print, tmp_path):
        output = tmp_path / "tree_space_kf.png"
        testargs = [
            "phykit", "tree_space",
            "-t", "tests/sample_files/gene_trees_simple.nwk",
            "-o", str(output),
            "--metric", "kf",
            "--seed", "42",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output_text = "\n".join(str(call) for call in mocked_print.call_args_list)
        assert "Metric: KF" in output_text
        assert output.exists()

    @patch("builtins.print")
    def test_tsne_method(self, mocked_print, tmp_path):
        output = tmp_path / "tree_space_tsne.png"
        testargs = [
            "phykit", "tree_space",
            "-t", "tests/sample_files/gene_trees_simple.nwk",
            "-o", str(output),
            "--method", "tsne",
            "--seed", "42",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output_text = "\n".join(str(call) for call in mocked_print.call_args_list)
        assert "Method: TSNE" in output_text
        assert output.exists()

    @patch("builtins.print")
    def test_species_tree(self, mocked_print, tmp_path):
        output = tmp_path / "tree_space_sp.png"
        testargs = [
            "phykit", "tree_space",
            "-t", "tests/sample_files/gene_trees_simple.nwk",
            "-o", str(output),
            "--species-tree", "tests/sample_files/tree_simple.tre",
            "--seed", "42",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output_text = "\n".join(str(call) for call in mocked_print.call_args_list)
        assert "Species tree: Cluster" in output_text
        assert output.exists()

    @patch("builtins.print")
    def test_json_output(self, mocked_print, tmp_path):
        output = tmp_path / "tree_space.png"
        testargs = [
            "phykit", "tree_space",
            "-t", "tests/sample_files/gene_trees_simple.nwk",
            "-o", str(output),
            "--seed", "42",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        json_str = mocked_print.call_args.args[0]
        payload = json.loads(json_str)
        assert payload["method"] == "mds"
        assert payload["metric"] == "rf"
        assert payload["n_trees"] == 10
        assert payload["n_taxa"] == 8
        assert "clusters" in payload
        assert "coordinates" in payload

    @patch("builtins.print")
    def test_seed_reproducible(self, mocked_print, tmp_path):
        # Run twice with same seed
        results = []
        for i in range(2):
            output = tmp_path / f"tree_space_{i}.png"
            testargs = [
                "phykit", "tree_space",
                "-t", "tests/sample_files/gene_trees_simple.nwk",
                "-o", str(output),
                "--seed", "42",
                "--json",
            ]
            with patch.object(sys, "argv", testargs):
                Phykit()

            json_str = mocked_print.call_args.args[0]
            payload = json.loads(json_str)
            results.append(payload["coordinates"])

        # Coordinates should be identical across runs
        assert results[0] == results[1]
