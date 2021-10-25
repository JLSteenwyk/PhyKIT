import pytest
import sys
from math import isclose
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestRootTree(object):
    @patch("builtins.print")
    def test_root_tree(self, mocked_print):
        testargs = [
            "phykit",
            "root_tree",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-r",
            f"{here.parent.parent.parent}/sample_files/tree_simple.outgroup.txt",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/tree_simple.tre.rooted", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent.parent}/sample_files/tree_simple.tre.rooted", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

    @patch("builtins.print")
    def test_root_tree_alias0(self, mocked_print):
        testargs = [
            "phykit",
            "root",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-r",
            f"{here.parent.parent.parent}/sample_files/tree_simple.outgroup.txt",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/tree_simple.tre.rooted", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent.parent}/sample_files/tree_simple.tre.rooted", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

    @patch("builtins.print")
    def test_root_tree_alias1(self, mocked_print):
        testargs = [
            "phykit",
            "rt",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-r",
            f"{here.parent.parent.parent}/sample_files/tree_simple.outgroup.txt",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/tree_simple.tre.rooted", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent.parent}/sample_files/tree_simple.tre.rooted", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

    @patch("builtins.print")
    def test_root_tree_incorrect_tree_path(self, mocked_print):
        testargs = [
            "phykit",
            "root_tree",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tr",
            "-r",
            f"{here.parent.parent.parent}/sample_files/tree_simple.outgroup.txt",
        ]

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_root_tree_incorrect_root_path(self, mocked_print):
        testargs = [
            "phykit",
            "root_tree",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-r",
            f"{here.parent.parent.parent}/sample_files/tree_simple.outgroup.tx",
        ]

        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type == SystemExit
        mocked_print.assert_has_calls([
            call("Please check file name and pathing"),
        ])

    @patch("builtins.print")
    def test_root_tree_custom_output(self, mocked_print):
        testargs = [
            "phykit",
            "root_tree",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-r",
            f"{here.parent.parent.parent}/sample_files/tree_simple.outgroup.txt",
            "-o",
            f"{here.parent.parent.parent}/sample_files/tree_simple_rooted_custom_out.tre"
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/tree_simple.tre.rooted", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent.parent}/sample_files/tree_simple_rooted_custom_out.tre", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content