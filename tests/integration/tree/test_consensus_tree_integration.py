import json
import sys

import pytest
from mock import patch

from phykit.phykit import Phykit


@pytest.mark.integration
class TestConsensusTreeIntegration:
    @patch("builtins.print")
    def test_consensus_tree_json(self, mocked_print, tmp_path):
        trees = tmp_path / "trees.nwk"
        trees.write_text("((A,B),(C,D));\n((A,B),(C,D));\n")

        testargs = ["phykit", "consensus_tree", "-t", str(trees), "--json"]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["method"] == "majority"
        assert payload["taxa_in_consensus"] == 4

    @patch("builtins.print")
    def test_consensus_tree_alias(self, mocked_print, tmp_path):
        trees = tmp_path / "trees.nwk"
        trees.write_text("((A,B),(C,D));\n((A,B),(C,D));\n")

        testargs = ["phykit", "consensus", "-t", str(trees)]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert "A" in mocked_print.call_args.args[0]

    @patch("builtins.print")
    def test_consensus_tree_missing_taxa_shared(self, mocked_print, tmp_path):
        trees = tmp_path / "trees_missing.nwk"
        trees.write_text("((A,B),(C,D));\n((A,B),(C,E));\n")

        testargs = [
            "phykit",
            "ctree",
            "-t",
            str(trees),
            "--missing-taxa",
            "shared",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["pruned_to_shared_taxa"] is True
        assert payload["taxa_in_consensus"] == 3
        assert "D" not in payload["tree_newick"]
        assert "E" not in payload["tree_newick"]
