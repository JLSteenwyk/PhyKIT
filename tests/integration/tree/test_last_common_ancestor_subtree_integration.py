import pytest
import sys
from mock import patch
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestLastCommonAncestorSubtree(object):
    @patch("builtins.print")
    def test_lca_subtree_simple(self, mocked_print):
        testargs = [
            "phykit",
            "last_common_ancestor_subtree",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_lca_subtree_list.txt",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/tree_simple.tre.subtree.tre", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent.parent}/sample_files/tree_simple.tre.subtree.tre", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

    @patch("builtins.print")
    def test_lca_subtree_simple_alias(self, mocked_print):
        testargs = [
            "phykit",
            "lca_subtree",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_lca_subtree_list.txt",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/tree_simple.tre.subtree.tre", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent.parent}/sample_files/tree_simple.tre.subtree.tre", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

    @patch("builtins.print")
    def test_lca_subtree_incorrect_tree_path(self, mocked_print):
        testargs = [
            "phykit",
            "lca_subtree",
            f"{here.parent.parent.parent}/sample_files/tree_simple.t",
            f"{here.parent.parent.parent}/sample_files/tree_simple_lca_subtree_list.txt",
        ]

        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_lca_subtree_incorrect_list_path(self, mocked_print):
        testargs = [
            "phykit",
            "lca_subtree",
            f"{here.parent.parent.parent}/sample_files/tree_simple.t",
            f"{here.parent.parent.parent}/sample_files/tree_simple_lca_subtree_list.tx",
        ]

        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2
        # mocked_print.assert_has_calls([
        #     call("Please check file name and pathing"),
        # ])

    @patch("builtins.print")
    def test_rename_tree_tips_custom_output(self, mocked_print):
        testargs = [
            "phykit",
            "lca_subtree",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_lca_subtree_list.txt",
            '-o',
            f"{here.parent.parent.parent}/sample_files/tree_simple_lca_subtree.tre"

        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/tree_simple_lca_subtree.tre", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent.parent}/sample_files/tree_simple_lca_subtree.tre", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content
