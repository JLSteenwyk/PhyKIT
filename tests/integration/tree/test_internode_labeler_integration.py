import pytest
import sys
import json
from mock import patch
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestInternodeLabeler(object):
    @patch("builtins.print")
    def test_internode_labeler(self, mocked_print):
        testargs = [
            "phykit",
            "internode_labeler",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/tree_simple.tre.internode_labels.tre", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent.parent}/sample_files/tree_simple.tre.internode_labels.tre", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

    @patch("builtins.print")
    def test_internode_labeler_alias(self, mocked_print):
        testargs = [
            "phykit",
            "il",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/tree_simple.tre.internode_labels.tre", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent.parent}/sample_files/tree_simple.tre.internode_labels.tre", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

    @patch("builtins.print")
    def test_internode_labeler_custom_out(self, mocked_print):
        testargs = [
            "phykit",
            "internode_labeler",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-o",
            f"{here.parent.parent.parent}/sample_files/tree_simple.custom_out_internode_labels.tre"
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/tree_simple.tre.internode_labels.tre", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent.parent}/sample_files/tree_simple.tre.internode_labels.tre", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

    @patch("builtins.print")
    def test_internode_labeler_incorrect_path(self, mocked_print):

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_internode_labeler_json(self, mocked_print):
        testargs = [
            "phykit",
            "internode_labeler",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["labeled_internodes"] == 6
        assert payload["output_file"].endswith(".internode_labels.tre")
