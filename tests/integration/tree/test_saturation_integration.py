import pytest
import sys
from math import isclose
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestSaturation(object):
    @patch("builtins.print")
    def test_saturation(self, mocked_print):
        expected_result = 0.8451
        testargs = [
            "phykit",
            "saturation",
            "-t",
            f"{here.parent.parent.parent}/sample_files/12_YPR191W_Anc_7.548_codon_aln.fasta.clipkit.treefile",
            "-a",
            f"{here.parent.parent.parent}/sample_files/12_YPR191W_Anc_7.548_codon_aln.fasta.clipkit",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_saturation_alias(self, mocked_print):
        expected_result = 0.8451
        testargs = [
            "phykit",
            "sat",
            "-t",
            f"{here.parent.parent.parent}/sample_files/12_YPR191W_Anc_7.548_codon_aln.fasta.clipkit.treefile",
            "-a",
            f"{here.parent.parent.parent}/sample_files/12_YPR191W_Anc_7.548_codon_aln.fasta.clipkit",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_saturation_incorrect_tree_path(self, mocked_print):
        testargs = [
            "phykit",
            "saturation",
            "-t",
            f"{here.parent.parent.parent}/sample_files/12_YPR191W_Anc_7.548_codon_aln.fasta.clipkit.treefi",
            "-a",
            f"{here.parent.parent.parent}/sample_files/12_YPR191W_Anc_7.548_codon_aln.fasta.clipkit",
        ]

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_saturation_incorrect_alignment_path(self, mocked_print):
        testargs = [
            "phykit",
            "saturation",
            "-t",
            f"{here.parent.parent.parent}/sample_files/12_YPR191W_Anc_7.548_codon_aln.fasta.clipkit.treefile",
            "-a",
            f"{here.parent.parent.parent}/sample_files/12_YPR191W_Anc_7.548_codon_aln.fasta.clipki",
        ]

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_saturation_verbose(self, mocked_print):
        expected_result = 0.8451
        testargs = [
            "phykit",
            "sat",
            "-t",
            f"{here.parent.parent.parent}/sample_files/12_YPR191W_Anc_7.548_codon_aln.fasta.clipkit.treefile",
            "-a",
            f"{here.parent.parent.parent}/sample_files/12_YPR191W_Anc_7.548_codon_aln.fasta.clipkit",
            "-v"
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call("Kpol-Kpha\t0.6136\t0.6176"),
            call("Kpol-Snag\t0.5654\t0.7482"),
            call("Kpol-Suva\t0.5948\t0.6945"),
            call("Kpol-Skud\t0.5906\t0.7341"),
            call("Kpol-Smik\t0.6157\t0.7734"),
            call("Kpol-Scer\t0.6105\t0.7634"),
            call("Kpol-Kbla\t0.599\t0.7287"),
            call("Kpol-Kafr\t0.5696\t0.7509"),
            call("Kpol-Sdai\t0.5916\t0.7297"),
            call("Kpol-Scas\t0.5822\t0.7958"),
            call("Kpol-Cgla\t0.5979\t0.653"),
            call("Kpha-Snag\t0.5393\t0.8254"),
            call("Kpha-Suva\t0.5424\t0.7717"),
            call("Kpha-Skud\t0.5257\t0.8113"),
            call("Kpha-Smik\t0.5466\t0.8505"),
            call("Kpha-Scer\t0.5445\t0.8406"),
            call("Kpha-Kbla\t0.5435\t0.9228"),
            call("Kpha-Kafr\t0.5581\t0.945"),
            call("Kpha-Sdai\t0.5424\t0.9238"),
            call("Kpha-Scas\t0.5382\t0.9899"),
            call("Kpha-Cgla\t0.5738\t0.8472"),
            call("Snag-Suva\t0.5518\t0.7281"),
            call("Snag-Skud\t0.5508\t0.7678"),
            call("Snag-Smik\t0.5539\t0.807"),
            call("Snag-Scer\t0.555\t0.797"),
            call("Snag-Kbla\t0.5215\t1.0535"),
            call("Snag-Kafr\t0.5539\t1.0757"),
            call("Snag-Sdai\t0.5225\t1.0544"),
            call("Snag-Scas\t0.4984\t1.1206"),
            call("Snag-Cgla\t0.5131\t0.9778"),
            call("Suva-Skud\t0.8105\t0.229"),
            call("Suva-Smik\t0.7969\t0.2682"),
            call("Suva-Scer\t0.7958\t0.2583"),
            call("Suva-Kbla\t0.5361\t0.9997"),
            call("Suva-Kafr\t0.5162\t1.0219"),
            call("Suva-Sdai\t0.5309\t1.0007"),
            call("Suva-Scas\t0.5162\t1.0668"),
            call("Suva-Cgla\t0.5707\t0.9241"),
            call("Skud-Smik\t0.8199\t0.2171"),
            call("Skud-Scer\t0.8325\t0.2071"),
            call("Skud-Kbla\t0.5183\t1.0394"),
            call("Skud-Kafr\t0.5288\t1.0616"),
            call("Skud-Sdai\t0.5089\t1.0403"),
            call("Skud-Scas\t0.5005\t1.1065"),
            call("Skud-Cgla\t0.5602\t0.9637"),
            call("Smik-Scer\t0.8314\t0.1998"),
            call("Smik-Kbla\t0.5382\t1.0786"),
            call("Smik-Kafr\t0.5225\t1.1008"),
            call("Smik-Sdai\t0.534\t1.0796"),
            call("Smik-Scas\t0.5257\t1.1457"),
            call("Smik-Cgla\t0.5529\t1.003"),
            call("Scer-Kbla\t0.5372\t1.0686"),
            call("Scer-Kafr\t0.5162\t1.0908"),
            call("Scer-Sdai\t0.5298\t1.0696"),
            call("Scer-Scas\t0.5236\t1.1357"),
            call("Scer-Cgla\t0.5508\t0.993"),
            call("Kbla-Kafr\t0.5686\t0.67"),
            call("Kbla-Sdai\t0.5613\t0.8756"),
            call("Kbla-Scas\t0.5393\t0.9418"),
            call("Kbla-Cgla\t0.5508\t0.799"),
            call("Kafr-Sdai\t0.534\t0.8978"),
            call("Kafr-Scas\t0.5152\t0.964"),
            call("Kafr-Cgla\t0.555\t0.8212"),
            call("Sdai-Scas\t0.6356\t0.5357"),
            call("Sdai-Cgla\t0.5675\t0.7045"),
            call("Scas-Cgla\t0.5518\t0.7706"),
        ]