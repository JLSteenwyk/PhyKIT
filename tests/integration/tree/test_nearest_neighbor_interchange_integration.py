import pytest
import sys
from mock import patch
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestNearestNeighborInterchange(object):
    @patch("builtins.print")
    def test_nni(self, mocked_print):
        testargs = [
            "phykit",
            "nearest_neighbor_interchange",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tre_rooted.tree",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/small_Aspergillus_tre_rooted.nni_moves.tre", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tre_rooted.tree.nnis", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

    @patch("builtins.print")
    def test_nni_custom_out(self, mocked_print):
        testargs = [
            "phykit",
            "nearest_neighbor_interchange",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tre_rooted.tree",
            "-o",
            f"{here.parent.parent.parent}/sample_files/nni_custom_out.tre"
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/small_Aspergillus_tre_rooted.nni_moves.tre", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent.parent}/sample_files/nni_custom_out.tre", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

    @patch("builtins.print")
    def test_nni_alias(self, mocked_print):
        testargs = [
            "phykit",
            "nni",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tre_rooted.tree",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/small_Aspergillus_tre_rooted.nni_moves.tre", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tre_rooted.tree.nnis", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

    @patch("builtins.print")
    def test_nni_wrong_path_tre(self, mocked_print):
        testargs = [
            "phykit",
            "nni",
            "bad_path",
        ]

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2
