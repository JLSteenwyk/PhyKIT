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
class TestPruneTree(object):
    @patch("builtins.print")
    def test_prune(self, mocked_print):
        testargs = [
            "phykit",
            "prune",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_prune.txt",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/tree_simple.tre.pruned", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent.parent}/sample_files/tree_simple.tre.pruned", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

    @patch("builtins.print")
    def test_prune_custom_out(self, mocked_print):
        testargs = [
            "phykit",
            "prune",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_prune.txt",
            "-o",
            f"{here.parent.parent.parent}/sample_files/tree_simple_pruned_custom_output.tre"
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/tree_simple_pruned_custom_output.tre", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent.parent}/sample_files/tree_simple_pruned_custom_output.tre", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

    @patch("builtins.print")
    def test_prune_wrong_path_tre(self, mocked_print):
        testargs = [
            "phykit",
            "prune",
            "/does/not/exist.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_prune.txt",
        ]
        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_prune_wrong_path_list(self, mocked_print):
        testargs = [
            "phykit",
            "prune",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "/does/not/exist.txt",
        ]
        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_prune_alias(self, mocked_print):
        testargs = [
            "phykit",
            "prune_tree",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_prune.txt",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/tree_simple.tre.pruned", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent.parent}/sample_files/tree_simple.tre.pruned", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

    @patch("builtins.print")
    def test_prune_keep(self, mocked_print):
        testargs = [
            "phykit",
            "prune",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_prune_keep.txt",
            "-k"
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/tree_simple.pruned_keep.tre", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent.parent}/sample_files/tree_simple.tre.pruned", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

    @patch("builtins.print")
    def test_prune_keep_long(self, mocked_print):
        testargs = [
            "phykit",
            "prune",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_prune_keep.txt",
            "--keep"
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/tree_simple.pruned_keep.tre", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent.parent}/sample_files/tree_simple.tre.pruned", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

    @patch("builtins.print")
    def test_prune_json(self, mocked_print):
        testargs = [
            "phykit",
            "prune",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_prune.txt",
            "--json",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["keep_input_taxa"] is False
        assert payload["pruned_count"] == 2
        assert payload["remaining_tips"] == 6

    @patch("builtins.print")
    def test_prune_keep_ignore_branch_labels(self, mocked_print):
        out_path = (
            f"{here.parent.parent.parent}/sample_files/"
            "tree_simple_labeled.tre.keep_stripped"
        )
        testargs = [
            "phykit",
            "prune",
            f"{here.parent.parent.parent}/sample_files/tree_simple_labeled.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_labeled_keep.txt",
            "-k",
            "--ignore-branch-labels",
            "-o",
            out_path,
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(out_path) as fh:
            tree = Phylo.read(StringIO(fh.read()), "newick")

        tip_names = sorted(t.name for t in tree.get_terminals())
        # Keep list (raccoon, bear, monkey, cat) matches by stripped name;
        # the {FG} markers must be preserved on the kept tips.
        assert tip_names == ["bear", "cat{FG}", "monkey{FG}", "raccoon"]

    @patch("builtins.print")
    def test_prune_ignore_branch_labels(self, mocked_print):
        out_path = (
            f"{here.parent.parent.parent}/sample_files/"
            "tree_simple_labeled.tre.prune_stripped"
        )
        testargs = [
            "phykit",
            "prune",
            f"{here.parent.parent.parent}/sample_files/tree_simple_labeled.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_labeled_prune.txt",
            "--ignore-branch-labels",
            "-o",
            out_path,
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(out_path) as fh:
            tree = Phylo.read(StringIO(fh.read()), "newick")

        tip_names = sorted(t.name for t in tree.get_terminals())
        # Prune list (monkey, cat) matches by stripped name; both labeled
        # tips should be removed even though the list lacks {FG}.
        assert tip_names == ["bear", "dog", "raccoon", "sea_lion", "seal", "weasel"]

    @patch("builtins.print")
    def test_prune_keep_ignore_branch_labels_json(self, mocked_print):
        out_path = (
            f"{here.parent.parent.parent}/sample_files/"
            "tree_simple_labeled.tre.keep_stripped_json"
        )
        testargs = [
            "phykit",
            "prune",
            f"{here.parent.parent.parent}/sample_files/tree_simple_labeled.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_labeled_keep.txt",
            "-k",
            "--ignore-branch-labels",
            "-o",
            out_path,
            "--json",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["keep_input_taxa"] is True
        assert payload["ignore_branch_labels"] is True
        assert payload["remaining_tips"] == 4
        # Pruned names retain their {FG} labels in the report.
        assert payload["taxa_pruned"] == sorted(
            ["sea_lion", "seal", "weasel", "dog"]
        )
