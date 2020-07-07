import pytest
import sys
from math import isclose
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestTree(object):
    @patch("builtins.print")
    def test_rename_tree_tips(self, mocked_print):
        testargs = [
            "phykit",
            "rename_tree_tips",
            f"{here.parent.parent}/sample_files/tree_simple.tre",
            "-i",
            f"{here.parent.parent}/sample_files/tree_simple_idmap.txt",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent}/expected/tree_simple.tre.renamed.tre", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent}/sample_files/tree_simple.tre.renamed.tre", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

    @patch("builtins.print")
    def test_rf_distance(self, mocked_print):
        expected_result = "8\t0.8"
        testargs = [
            "phykit",
            "rf_distance",
            f"{here.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent}/sample_files/tree_simple_other_topology.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]
    
    @patch("builtins.print")
    def test_saturation(self, mocked_print):
        expected_result = 0.8451
        testargs = [
            "phykit",
            "saturation",
            "-t",
            f"{here.parent.parent}/sample_files/12_YPR191W_Anc_7.548_codon_aln.fasta.clipkit.treefile",
            "-a",
            f"{here.parent.parent}/sample_files/12_YPR191W_Anc_7.548_codon_aln.fasta.clipkit",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_spurious_sequence(self, mocked_print):
        expected_result = "None"
        testargs = [
            "phykit",
            "spurious_sequence",
            f"{here.parent.parent}/sample_files/tree_simple.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_tip_labels(self, mocked_print):
        expected_result = "None"
        testargs = [
            "phykit",
            "tip_labels",
            f"{here.parent.parent}/sample_files/tree_simple.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call("raccoon"),
            call("bear"),
            call("sea_lion"),
            call("seal"),
            call("monkey"),
            call("cat"),
            call("weasel"),
            call("dog"),
        ]

    @patch("builtins.print")
    def test_total_tree_length(self, mocked_print):
        expected_result = 277.2772
        testargs = [
            "phykit",
            "total_tree_length",
            f"{here.parent.parent}/sample_files/tree_simple.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_treeness(self, mocked_print):
        expected_result = 0.1899
        testargs = [
            "phykit",
            "treeness",
            f"{here.parent.parent}/sample_files/Yeasts_2832_eMRC_reference_renamed.tree",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_treeness_over_rcv(self, mocked_print):
        expected_result = "0.35\t0.126\t0.36"
        testargs = [
            "phykit",
            "treeness_over_rcv",
            "-t",
            f"{here.parent.parent}/sample_files/tree_simple.tre",
            "-a",
            f"{here.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]