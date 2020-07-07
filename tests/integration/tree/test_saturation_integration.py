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