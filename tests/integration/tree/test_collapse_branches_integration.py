import pytest
import sys
from mock import patch
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestCollapseBranches(object):
    @patch("builtins.print")
    def test_collapse_branches0(self, mocked_print):
        testargs = [
            "phykit",
            "collapse_branches",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            "-s",
            "100",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/small_Aspergillus_tree.tre.collapsed_100.0.tre", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre.collapsed_100.0.tre", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

    @patch("builtins.print")
    def test_collapse_branches1(self, mocked_print):
        testargs = [
            "phykit",
            "collapse_branches",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            "-s",
            "80",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/small_Aspergillus_tree.tre.collapsed_80.0.tre", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre.collapsed_80.0.tre", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

    @patch("builtins.print")
    def test_collapse_branches2(self, mocked_print):
        testargs = [
            "phykit",
            "collapse_branches",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            "-s",
            "90",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/small_Aspergillus_tree.tre.collapsed_90.0.tre", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre.collapsed_90.0.tre", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

    @patch("builtins.print")
    def test_collapse_branches_incorrect_input(self, mocked_print):
        testargs = [
            "phykit",
            "collapse_branches",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tr",
            "-s",
            "90",
        ]

        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_collapse_branches_incorrect_args(self, mocked_print):
        testargs = [
            "phykit",
            "collapse_branches",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            "-s",
        ]

        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_collapse_branches_alias0(self, mocked_print):
        testargs = [
            "phykit",
            "collapse",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            "-s",
            "90",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/small_Aspergillus_tree.tre.collapsed_90.0.tre", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre.collapsed_90.0.tre", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

    @patch("builtins.print")
    def test_collapse_branches_alias1(self, mocked_print):
        testargs = [
            "phykit",
            "cb",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            "-s",
            "90",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/small_Aspergillus_tree.tre.collapsed_90.0.tre", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre.collapsed_90.0.tre", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content
