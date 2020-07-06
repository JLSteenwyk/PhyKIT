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
    def test_dna_threader0(self, mocked_print):
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

    @patch("builtins.print")
    def test_dna_threader1(self, mocked_print):
        expected_result_0 = dedent(
            """>1\nAAAGGG---"""
        )
        expected_result_1 = dedent(
            """>2\nAAATTTGGG"""  
        )
        expected_result_2 = dedent(
            """>3\nAAATTTGGG"""  
        )
        expected_result_3 = dedent(
            """>4\nAAATTTGGG"""  
        )
        testargs = [
            "phykit",
            "thread_dna",
            "-p",
            f"{here.parent.parent.parent}/sample_files/test_alignment.prot.faa",
            "-n",
            f"{here.parent.parent.parent}/sample_files/test.nucl.fna",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
            call(expected_result_2),
            call(expected_result_3)
        ]

    @patch("builtins.print")
    def test_dna_threader_alias0(self, mocked_print):
        expected_result_0 = dedent(
            """>1\nAAAGGG---"""
        )
        expected_result_1 = dedent(
            """>2\nAAATTTGGG"""  
        )
        expected_result_2 = dedent(
            """>3\nAAATTTGGG"""  
        )
        expected_result_3 = dedent(
            """>4\nAAATTTGGG"""  
        )
        testargs = [
            "phykit",
            "pal2nal",
            "-p",
            f"{here.parent.parent.parent}/sample_files/test_alignment.prot.faa",
            "-n",
            f"{here.parent.parent.parent}/sample_files/test.nucl.fna",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
            call(expected_result_2),
            call(expected_result_3)
        ]

    @patch("builtins.print")
    def test_dna_threader_alias1(self, mocked_print):
        expected_result_0 = dedent(
            """>1\nAAAGGG---"""
        )
        expected_result_1 = dedent(
            """>2\nAAATTTGGG"""  
        )
        expected_result_2 = dedent(
            """>3\nAAATTTGGG"""  
        )
        expected_result_3 = dedent(
            """>4\nAAATTTGGG"""  
        )
        testargs = [
            "phykit",
            "p2n",
            "-p",
            f"{here.parent.parent.parent}/sample_files/test_alignment.prot.faa",
            "-n",
            f"{here.parent.parent.parent}/sample_files/test.nucl.fna",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
            call(expected_result_2),
            call(expected_result_3)
        ]

    @patch("builtins.print")
    def test_dna_threader_incorrect_input_file(self, mocked_print):
        testargs = [
            "phykit",
            "thread_dna",
            "-p",
            f"{here.parent.parent.parent}/sample_files/test_alignment.prot.fa",
            "-n",
            f"{here.parent.parent.parent}/sample_files/test.nucl.fna",
        ]

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_dna_threader_alias_stop_codon_true(self, mocked_print):
        expected_result_0 = dedent(
            """>1\nAAAGGG---"""
        )
        expected_result_1 = dedent(
            """>2\nAAATTTGGG"""  
        )
        expected_result_2 = dedent(
            """>3\nAAATTTGGG"""  
        )
        expected_result_3 = dedent(
            """>4\nAAATTT---"""  
        )
        testargs = [
            "phykit",
            "p2n",
            "-p",
            f"{here.parent.parent.parent}/sample_files/test_alignment.prot.faa",
            "-n",
            f"{here.parent.parent.parent}/sample_files/test.nucl.fna",
            "-s"
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
            call(expected_result_2),
            call(expected_result_3)
        ]