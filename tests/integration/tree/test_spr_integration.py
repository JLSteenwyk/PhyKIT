import pytest
import sys
import json
from io import StringIO
from mock import patch
from pathlib import Path

from Bio import Phylo

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestSprIntegration(object):
    @patch("builtins.print")
    def test_basic_invocation(self, mocked_print):
        """End-to-end test: subtree_prune_regraft with -t and --subtree."""
        testargs = [
            "phykit",
            "subtree_prune_regraft",
            "-t",
            f"{here.parent.parent.parent}/sample_files/spr_test.tree",
            "--subtree",
            "D,E",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        # Should print Newick trees to stdout
        assert mocked_print.call_count >= 1
        # Each call should be a valid Newick string
        for call in mocked_print.call_args_list:
            newick_str = call.args[0]
            tree = Phylo.read(StringIO(newick_str), "newick")
            taxa = {t.name for t in tree.get_terminals()}
            assert taxa == {"A", "B", "C", "D", "E"}

    @patch("builtins.print")
    def test_alias_spr(self, mocked_print):
        """The 'spr' alias works."""
        testargs = [
            "phykit",
            "spr",
            "-t",
            f"{here.parent.parent.parent}/sample_files/spr_test.tree",
            "--subtree",
            "D,E",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.call_count >= 1

    @patch("builtins.print")
    def test_json_output(self, mocked_print):
        """JSON output has correct structure."""
        testargs = [
            "phykit",
            "subtree_prune_regraft",
            "-t",
            f"{here.parent.parent.parent}/sample_files/spr_test.tree",
            "--subtree",
            "D,E",
            "--json",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        # The JSON should be printed in a single call
        assert mocked_print.call_count >= 1
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["subtree_taxa"] == ["D", "E"]
        assert payload["n_spr_trees"] > 0
        assert len(payload["trees"]) == payload["n_spr_trees"]
        for entry in payload["trees"]:
            assert "description" in entry
            assert "newick" in entry

    @patch("builtins.print")
    def test_output_file(self, mocked_print, tmp_path):
        """The -o flag creates a file with multiple trees."""
        out_path = str(tmp_path / "spr_trees.nwk")
        testargs = [
            "phykit",
            "subtree_prune_regraft",
            "-t",
            f"{here.parent.parent.parent}/sample_files/spr_test.tree",
            "--subtree",
            "D,E",
            "-o",
            out_path,
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(out_path) as f:
            content = f.read()
        lines = [l for l in content.strip().split("\n") if l.strip()]
        assert len(lines) >= 1
        # Each line should be valid Newick
        for line in lines:
            tree = Phylo.read(StringIO(line), "newick")
            taxa = {t.name for t in tree.get_terminals()}
            assert taxa == {"A", "B", "C", "D", "E"}

    @patch("builtins.print")
    def test_single_taxon_subtree(self, mocked_print):
        """SPR with a single taxon subtree works."""
        testargs = [
            "phykit",
            "spr",
            "-t",
            f"{here.parent.parent.parent}/sample_files/spr_test.tree",
            "--subtree",
            "C",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.call_count >= 1

    @patch("builtins.print")
    def test_missing_taxon_raises(self, mocked_print):
        """Nonexistent taxon raises error with exit code 2."""
        testargs = [
            "phykit",
            "spr",
            "-t",
            f"{here.parent.parent.parent}/sample_files/spr_test.tree",
            "--subtree",
            "NONEXISTENT",
        ]

        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_bad_tree_path(self, mocked_print):
        """Nonexistent tree file raises error with exit code 2."""
        testargs = [
            "phykit",
            "spr",
            "-t",
            "no_such_file.tree",
            "--subtree",
            "A",
        ]

        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 2
