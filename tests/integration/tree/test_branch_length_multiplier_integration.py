import pytest
import sys
from math import isclose
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestBranchLengthMultiplier(object):
    @patch("builtins.print")
    def test_branch_length_multiplier_custom_output(self, mocked_print):
        testargs = [
            "phykit",
            "branch_length_multiplier",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-f",
            "5",
            "-o",
            "./tests/sample_files/tree_simple_blm_5.tre"
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/tree_simple_blm_5.tre", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent.parent}/sample_files/tree_simple_blm_5.tre", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

    @patch("builtins.print")
    def test_branch_length_multiplier_default_output(self, mocked_print):
        testargs = [
            "phykit",
            "branch_length_multiplier",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-f",
            "2",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/tree_simple.tre.factor_2.0.tre", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent.parent}/sample_files/tree_simple.tre.factor_2.0.tre", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

    @patch("builtins.print")
    def test_branch_length_multiplier_alias(self, mocked_print):
        testargs = [
            "phykit",
            "blm",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-f",
            "2",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/tree_simple.tre.factor_2.0.tre", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent.parent}/sample_files/tree_simple.tre.factor_2.0.tre", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

    @patch("builtins.print")
    def test_branch_length_multiplier_incorrect_input(self, mocked_print):
        testargs = [
            "phykit",
            "blm",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tr",
            "-f",
            "2",
        ]
        
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_branch_length_multiplier_incorrect_factor(self, mocked_print):
        testargs = [
            "phykit",
            "blm",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tr",
            "-f",
        ]
        
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2