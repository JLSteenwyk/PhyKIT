import pytest
import sys
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestDNAThreader(object):
    @patch("builtins.print")
    def test_dna_threader(self, mocked_print):
        expected_result_0 = dedent(
            """>200_S38\natggctgacatcctcacgcagctccagacttgcctggatcagcttgcaacacaattctacgcaacacttggttatctcacaacataccacgacaatgcccccacaacaccaccacca------aatgtccccgacgcagcaccagccctagcaaagatcaccaagaactcatcatcaccgccagtcccagcagccatcgcaaataaagtggggggtgcagctgctgttgcgggcaatgca---tcacccccacaggcgcct---------cct---------------------------------------------caacaacccggagct---------gcgcca---gggagagcagtagaaggtgaagatcccaaccttcctcccgcgccagactcgcccagcacgtttgcaagccggcagcgggagcttgcgcgcgatctcattatcaaagaacagcagatcgagtaccttatctccgtgcttcccgggattggcgcctctgaggctgaacaagaaaccagaatccaggacctggagaccgagcttagagacgtcgagaaggagcgcgctgcgaaagtgcgggagttgaaaaagttgaggactcggttggaggatgttcttggcgctgtcgctgtgggtatccacggggatggttactctcaaaac---------"""
        )
        expected_result_1 = dedent(
            """>203_S40\natggctgacatcctcacgcagctccagacttgcctggatcagcttgcaacacaattctacgcaacacttggttatctcacaacataccacgacaatgcccccacaacaccaccacca------aatgtccccgacgcagcaccagccctagcaaagatcaccaagaactcatcatcaccaccagtcccagcagccatcgcaaataaagtggggggtgcagctgctgttgcgggcaatgca---tcacccccacaggcgcct---------cct---------------------------------------------caacaacccggagct---------gcgcca---gggagagcagtagaaggtgaagatcccaaccttcctcccgcgccagactcgcccagcacgtttgcaagccggcagcgggagcttgcgcgcgatctcattatcaaagaacagcagatcgagtaccttatctccgtgcttcccgggattggcgcctctgaggctgaacaagaaaccagaatccaggacctggagaccgagcttagagacgtcgagaaggagcgcgctgcgaaagtgcgggagttgaaaaagttgaggactcggttggaggatgttcttggcgctgtcgctgtgggtatccacggggatggttactctcaaaac---------"""
        )

        testargs = [
            "phykit",
            "thread_dna",
            "-p",
            f"{here.parent.parent.parent}/sample_files/EOG091N44MS.fa.mafft",
            "-n",
            f"{here.parent.parent.parent}/sample_files/EOG091N44MS.fa",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
        ]