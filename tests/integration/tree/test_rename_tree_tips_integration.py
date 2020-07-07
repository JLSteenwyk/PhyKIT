import pytest
import sys
from math import isclose
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestRenameTreeTips(object):
    @patch("builtins.print")
    def test_rename_tree_tips(self, mocked_print):
        testargs = [
            "phykit",
            "rename_tree_tips",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-i",
            f"{here.parent.parent.parent}/sample_files/tree_simple_idmap.txt",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/tree_simple.tre.renamed.tre", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent.parent}/sample_files/tree_simple.tre.renamed.tre", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

    @patch("builtins.print")
    def test_rename_tree_tips_alias0(self, mocked_print):
        testargs = [
            "phykit",
            "rename_tree",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-i",
            f"{here.parent.parent.parent}/sample_files/tree_simple_idmap.txt",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/tree_simple.tre.renamed.tre", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent.parent}/sample_files/tree_simple.tre.renamed.tre", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

    @patch("builtins.print")
    def test_rename_tree_tips_alias1(self, mocked_print):
        testargs = [
            "phykit",
            "rename_tips",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-i",
            f"{here.parent.parent.parent}/sample_files/tree_simple_idmap.txt",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/tree_simple.tre.renamed.tre", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent.parent}/sample_files/tree_simple.tre.renamed.tre", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

    @patch("builtins.print")
    def test_rename_tree_tips_incorrect_tree_path(self, mocked_print):
        testargs = [
            "phykit",
            "rename_tips",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tr",
            "-i",
            f"{here.parent.parent.parent}/sample_files/tree_simple_idmap.txt",
        ]

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_rename_tree_tips_incorrect_idmap_path(self, mocked_print):
        testargs = [
            "phykit",
            "rename_tips",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-i",
            f"{here.parent.parent.parent}/sample_files/tree_simple_idmap.tx",
        ]

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2