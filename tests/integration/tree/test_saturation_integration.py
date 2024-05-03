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
        expected_result = 0.4919
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
        expected_result = 0.4919
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
            f"{here.parent.parent.parent}/sample_files/does_not_exist",
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
            f"{here.parent.parent.parent}/sample_files/does_not_exist",
        ]

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_saturation_verbose(self, mocked_print):
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
            call("Kpol-Kpha\t0.3864\t0.6176"),
            call("Kpol-Snag\t0.4346\t0.7482"),
            call("Kpol-Suva\t0.4052\t0.6945"),
            call("Kpol-Skud\t0.4094\t0.7341"),
            call("Kpol-Smik\t0.3843\t0.7734"),
            call("Kpol-Scer\t0.3895\t0.7634"),
            call("Kpol-Kbla\t0.401\t0.7287"),
            call("Kpol-Kafr\t0.4304\t0.7509"),
            call("Kpol-Sdai\t0.4084\t0.7297"),
            call("Kpol-Scas\t0.4178\t0.7958"),
            call("Kpol-Cgla\t0.4021\t0.653"),
            call("Kpha-Snag\t0.4607\t0.8254"),
            call("Kpha-Suva\t0.4576\t0.7717"),
            call("Kpha-Skud\t0.4743\t0.8113"),
            call("Kpha-Smik\t0.4534\t0.8505"),
            call("Kpha-Scer\t0.4555\t0.8406"),
            call("Kpha-Kbla\t0.4565\t0.9228"),
            call("Kpha-Kafr\t0.4419\t0.945"),
            call("Kpha-Sdai\t0.4576\t0.9238"),
            call("Kpha-Scas\t0.4618\t0.9899"),
            call("Kpha-Cgla\t0.4262\t0.8472"),
            call("Snag-Suva\t0.4482\t0.7281"),
            call("Snag-Skud\t0.4492\t0.7678"),
            call("Snag-Smik\t0.4461\t0.807"),
            call("Snag-Scer\t0.445\t0.797"),
            call("Snag-Kbla\t0.4785\t1.0535"),
            call("Snag-Kafr\t0.4461\t1.0757"),
            call("Snag-Sdai\t0.4775\t1.0544"),
            call("Snag-Scas\t0.5016\t1.1206"),
            call("Snag-Cgla\t0.4869\t0.9778"),
            call("Suva-Skud\t0.1895\t0.229"),
            call("Suva-Smik\t0.2031\t0.2682"),
            call("Suva-Scer\t0.2042\t0.2583"),
            call("Suva-Kbla\t0.4639\t0.9997"),
            call("Suva-Kafr\t0.4838\t1.0219"),
            call("Suva-Sdai\t0.4691\t1.0007"),
            call("Suva-Scas\t0.4838\t1.0668"),
            call("Suva-Cgla\t0.4293\t0.9241"),
            call("Skud-Smik\t0.1801\t0.2171"),
            call("Skud-Scer\t0.1675\t0.2071"),
            call("Skud-Kbla\t0.4817\t1.0394"),
            call("Skud-Kafr\t0.4712\t1.0616"),
            call("Skud-Sdai\t0.4911\t1.0403"),
            call("Skud-Scas\t0.4995\t1.1065"),
            call("Skud-Cgla\t0.4398\t0.9637"),
            call("Smik-Scer\t0.1686\t0.1998"),
            call("Smik-Kbla\t0.4618\t1.0786"),
            call("Smik-Kafr\t0.4775\t1.1008"),
            call("Smik-Sdai\t0.466\t1.0796"),
            call("Smik-Scas\t0.4743\t1.1457"),
            call("Smik-Cgla\t0.4471\t1.003"),
            call("Scer-Kbla\t0.4628\t1.0686"),
            call("Scer-Kafr\t0.4838\t1.0908"),
            call("Scer-Sdai\t0.4702\t1.0696"),
            call("Scer-Scas\t0.4764\t1.1357"),
            call("Scer-Cgla\t0.4492\t0.993"),
            call("Kbla-Kafr\t0.4314\t0.67"),
            call("Kbla-Sdai\t0.4387\t0.8756"),
            call("Kbla-Scas\t0.4607\t0.9418"),
            call("Kbla-Cgla\t0.4492\t0.799"),
            call("Kafr-Sdai\t0.466\t0.8978"),
            call("Kafr-Scas\t0.4848\t0.964"),
            call("Kafr-Cgla\t0.445\t0.8212"),
            call("Sdai-Scas\t0.3644\t0.5357"),
            call("Sdai-Cgla\t0.4325\t0.7045"),
            call("Scas-Cgla\t0.4482\t0.7706"),
        ]
