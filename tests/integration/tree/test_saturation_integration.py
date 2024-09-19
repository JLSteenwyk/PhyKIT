from mock import patch, call
from pathlib import Path
import pytest
import sys

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestSaturation(object):
    @patch("builtins.print")
    def test_saturation(self, mocked_print):
        expected_result = "0.4919\t0.5081"
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
        expected_result = "0.4919\t0.5081"
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
            call("Kpol\tKpha\t0.3864\t0.6176"),
            call("Kpol\tSnag\t0.4346\t0.7482"),
            call("Kpol\tSuva\t0.4052\t0.6945"),
            call("Kpol\tSkud\t0.4094\t0.7341"),
            call("Kpol\tSmik\t0.3843\t0.7734"),
            call("Kpol\tScer\t0.3895\t0.7634"),
            call("Kpol\tKbla\t0.401\t0.7287"),
            call("Kpol\tKafr\t0.4304\t0.7509"),
            call("Kpol\tSdai\t0.4084\t0.7297"),
            call("Kpol\tScas\t0.4178\t0.7958"),
            call("Kpol\tCgla\t0.4021\t0.653"),
            call("Kpha\tSnag\t0.4607\t0.8254"),
            call("Kpha\tSuva\t0.4576\t0.7717"),
            call("Kpha\tSkud\t0.4743\t0.8113"),
            call("Kpha\tSmik\t0.4534\t0.8505"),
            call("Kpha\tScer\t0.4555\t0.8406"),
            call("Kpha\tKbla\t0.4565\t0.9228"),
            call("Kpha\tKafr\t0.4419\t0.945"),
            call("Kpha\tSdai\t0.4576\t0.9238"),
            call("Kpha\tScas\t0.4618\t0.9899"),
            call("Kpha\tCgla\t0.4262\t0.8472"),
            call("Snag\tSuva\t0.4482\t0.7281"),
            call("Snag\tSkud\t0.4492\t0.7678"),
            call("Snag\tSmik\t0.4461\t0.807"),
            call("Snag\tScer\t0.445\t0.797"),
            call("Snag\tKbla\t0.4785\t1.0535"),
            call("Snag\tKafr\t0.4461\t1.0757"),
            call("Snag\tSdai\t0.4775\t1.0544"),
            call("Snag\tScas\t0.5016\t1.1206"),
            call("Snag\tCgla\t0.4869\t0.9778"),
            call("Suva\tSkud\t0.1895\t0.229"),
            call("Suva\tSmik\t0.2031\t0.2682"),
            call("Suva\tScer\t0.2042\t0.2583"),
            call("Suva\tKbla\t0.4639\t0.9997"),
            call("Suva\tKafr\t0.4838\t1.0219"),
            call("Suva\tSdai\t0.4691\t1.0007"),
            call("Suva\tScas\t0.4838\t1.0668"),
            call("Suva\tCgla\t0.4293\t0.9241"),
            call("Skud\tSmik\t0.1801\t0.2171"),
            call("Skud\tScer\t0.1675\t0.2071"),
            call("Skud\tKbla\t0.4817\t1.0394"),
            call("Skud\tKafr\t0.4712\t1.0616"),
            call("Skud\tSdai\t0.4911\t1.0403"),
            call("Skud\tScas\t0.4995\t1.1065"),
            call("Skud\tCgla\t0.4398\t0.9637"),
            call("Smik\tScer\t0.1686\t0.1998"),
            call("Smik\tKbla\t0.4618\t1.0786"),
            call("Smik\tKafr\t0.4775\t1.1008"),
            call("Smik\tSdai\t0.466\t1.0796"),
            call("Smik\tScas\t0.4743\t1.1457"),
            call("Smik\tCgla\t0.4471\t1.003"),
            call("Scer\tKbla\t0.4628\t1.0686"),
            call("Scer\tKafr\t0.4838\t1.0908"),
            call("Scer\tSdai\t0.4702\t1.0696"),
            call("Scer\tScas\t0.4764\t1.1357"),
            call("Scer\tCgla\t0.4492\t0.993"),
            call("Kbla\tKafr\t0.4314\t0.67"),
            call("Kbla\tSdai\t0.4387\t0.8756"),
            call("Kbla\tScas\t0.4607\t0.9418"),
            call("Kbla\tCgla\t0.4492\t0.799"),
            call("Kafr\tSdai\t0.466\t0.8978"),
            call("Kafr\tScas\t0.4848\t0.964"),
            call("Kafr\tCgla\t0.445\t0.8212"),
            call("Sdai\tScas\t0.3644\t0.5357"),
            call("Sdai\tCgla\t0.4325\t0.7045"),
            call("Scas\tCgla\t0.4482\t0.7706"),
        ]
