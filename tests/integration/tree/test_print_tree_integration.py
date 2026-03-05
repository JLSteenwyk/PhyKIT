from pathlib import Path
import pytest
from mock import patch
import sys
import json

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestPrintTree(object):
    @patch("builtins.print")
    def test_print_tree0(self, mocked_print):
        testargs = [
            "phykit",
            "print_tree",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()


    @patch("builtins.print")
    def test_print_tree1(self, mocked_print):
        testargs = [
            "phykit",
            "print_tree",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()


    @patch("builtins.print")
    def test_print_tree_wrong_input(self, mocked_print):

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_print_tree_alias0(self, mocked_print):


        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 2


    @patch("builtins.print")
    def test_print_tree_alias1(self, mocked_print):


        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_print_tree_remove_branch_lengths_short(self, mocked_print):


        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_print_tree_remove_branch_lengths_long(self, mocked_print):


        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_print_tree_json(self, mocked_print):
        testargs = [
            "phykit",
            "print_tree",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["remove_branch_lengths"] is False
        assert "raccoon" in payload["tree_newick"]
        assert ":" in payload["tree_newick"]

    @patch("builtins.print")
    def test_print_tree_json_remove_branch_lengths(self, mocked_print):
        testargs = [
            "phykit",
            "print_tree",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "--remove",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["remove_branch_lengths"] is True
        assert "raccoon" in payload["tree_newick"]
        assert ":0.00000" in payload["tree_newick"]
