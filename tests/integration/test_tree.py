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
    def test_collapse_branches(self, mocked_print):
        testargs = [
            "phykit",
            "collapse_branches",
            f"{here.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            "-s",
            "100",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent}/expected/small_Aspergillus_tree.tre.collapsed_100.0.tre", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent}/sample_files/small_Aspergillus_tree.tre.collapsed_100.0.tre", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

    @patch("builtins.print")
    def test_covarying_evolutionary_rates(self, mocked_print):
        expected_result = "0.6769\t0.065228"
        testargs = [
            "phykit",
            "covarying_evolutionary_rates",
            f"{here.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent}/sample_files/tree_simple_1.tre",
            "-r",
            f"{here.parent.parent}/sample_files/tree_simple_2.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]
    
    @patch("builtins.print")
    def test_dvmc(self, mocked_print):
        expected_result = 42.8016
        testargs = [
            "phykit",
            "dvmc",
            "-t",
            f"{here.parent.parent}/sample_files/tree_simple.tre",
            "-r",
            f"{here.parent.parent}/sample_files/tree_simple.outgroup.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_internal_branch_stats(self, mocked_print):
        testargs = [
            "phykit",
            "internal_branch_stats",
            f"{here.parent.parent}/sample_files/tree_simple.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call("mean: 6.9872"),
            call("median: 3.8738"),
            call("25th percentile: 2.0946"),
            call("75th percentile: 7.5297"),
            call("minimum: 0.846"),
            call("maximum: 20.592"),
            call("standard deviation: 8.0114"),
            call("variance: 64.1826")
        ]

    @patch("builtins.print")
    def test_internode_labeler(self, mocked_print):
        testargs = [
            "phykit",
            "internode_labeler",
            f"{here.parent.parent}/sample_files/tree_simple.tre",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent}/expected/tree_simple.tre.internode_labels.tre", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent}/sample_files/tree_simple.tre.internode_labels.tre", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

    @patch("builtins.print")
    def test_lb_score(self, mocked_print):
        testargs = [
            "phykit",
            "lb_score",
            f"{here.parent.parent}/sample_files/tree_simple.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        expected = [
            dict(label="mean", value=-12.5),
            dict(label="median", value=-27.806),
            dict(label="25th percentile", value=-31.0492),
            dict(label="75th percentile", value=-12.9039),
            dict(label="minimum", value=-39.2834),
            dict(label="maximum", value=65.6709),
            dict(label="standard deviation", value=35.2669),
            dict(label="variance", value=1243.7527),
        ]

        for print_call, expected_call in zip(mocked_print.call_args_list, expected):
            print_call_args, _ = print_call
            [print_label, print_value] = print_call_args[0].split(': ')
            assert print_label == expected_call["label"]
            assert isclose(float(print_value), expected_call["value"], rel_tol = 0.0001)

    @patch("builtins.print")
    def test_patristic_distances(self, mocked_print):
        testargs = [
            "phykit",
            "patristic_distances",
            f"{here.parent.parent}/sample_files/tree_simple.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call("mean: 76.1974"),
            call("median: 49.5888"),
            call("25th percentile: 40.5054"),
            call("75th percentile: 108.1385"),
            call("minimum: 24.0"),
            call("maximum: 152.8813"),
            call("standard deviation: 45.4698"),
            call("variance: 2067.502")
        ]

    @patch("builtins.print")
    def test_polytomy_test(self, mocked_print):
        testargs = [
            "phykit",
            "polytomy_test",
            "-t",
            f"{here.parent.parent}/sample_files/test_trees.txt",
            "-g",
            f"{here.parent.parent}/sample_files/test_trees_groups.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call("Gene Support Frequency Results"),
            call("=============================="),
            call("chi-squared: 20.0"),
            call("p-value: 4.5e-05"),
            call("total genes: 10"),
            call("0-1: 10"),
            call("0-2: 0"),
            call("1-2: 0"),
            call("\nTriplet Results"),
            call("==============="),
            call("chi-squared: 400.0"),
            call("p-value: 0.0"),
            call("total triplets: 200"),
            call("0-1: 200"),
            call("0-2: 0"),
            call("1-2: 0")
        ]

    @patch("builtins.print")
    def test_prune(self, mocked_print):
        testargs = [
            "phykit",
            "prune",
            f"{here.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent}/sample_files/tree_simple_prune.txt",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent}/expected/tree_simple.tre.pruned", "r") as expected_tree:
            expected_tree_content = expected_tree.read()

        with open(f"{here.parent.parent}/sample_files/tree_simple.tre.pruned", "r") as out_tree:
            out_tree_content = out_tree.read()

        assert expected_tree_content == out_tree_content

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