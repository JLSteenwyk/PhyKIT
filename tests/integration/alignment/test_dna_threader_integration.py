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
            """>200_S38"""
        )
        expected_result_1 = dedent(
            """atggctgacatcctcacgcagctccagacttgcctggatcagcttgcaacacaattctac"""
        )
        expected_result_2 = dedent(
            """gcaacacttggttatctcacaacataccacgacaatgcccccacaacaccaccacca---"""
        )
        expected_result_3 = dedent(
            """---aatgtccccgacgcagcaccagccctagcaaagatcaccaagaactcatcatcaccg"""
        )
        expected_result_4 = dedent(
            """ccagtcccagcagccatcgcaaataaagtggggggtgcagctgctgttgcgggcaatgca"""
        )
        expected_result_5 = dedent(
            """>203_S40"""
        )
        expected_result_6 = dedent(
            """atggctgacatcctcacgcagctccagacttgcctggatcagcttgcaacacaattctac"""
        )
        expected_result_7 = dedent(
            """gcaacacttggttatctcacaacataccacgacaatgcccccacaacaccaccacca---"""
        )
        expected_result_8 = dedent(
            """---aatgtccccgacgcagcaccagccctagcaaagatcaccaagaactcatcatcacca"""
        )
        expected_result_9 = dedent(
            """ccagtcccagcagccatcgcaaataaagtggggggtgcagctgctgttgcgggcaatgca"""
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
            call(expected_result_2),
            call(expected_result_3),
            call(expected_result_4),
            call(expected_result_5),
            call(expected_result_6),
            call(expected_result_7),
            call(expected_result_8),
            call(expected_result_9),
        ]

    @patch("builtins.print")
    def test_dna_threader1(self, mocked_print):
        expected_result_0 = dedent(
            """>1"""
        )
        expected_result_1 = dedent(
            """AAAGGG---"""
        )
        expected_result_2 = dedent(
            """>2"""  
        )
        expected_result_3 = dedent(
            """AAATTTGGG"""
        )
        expected_result_4 = dedent(
            """>3"""  
        )
        expected_result_5 = dedent(
            """AAATTTGGG"""
        )
        expected_result_6 = dedent(
            """>4"""  
        )
        expected_result_7 = dedent(
            """AAATTTGGG"""
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
            call(expected_result_3),
            call(expected_result_4),
            call(expected_result_5),
            call(expected_result_6),
            call(expected_result_7),
        ]

    @patch("builtins.print")
    def test_dna_threader_alias0(self, mocked_print):
        expected_result_0 = dedent(
            """>1"""
        )
        expected_result_1 = dedent(
            """AAAGGG---"""
        )
        expected_result_2 = dedent(
            """>2"""  
        )
        expected_result_3 = dedent(
            """AAATTTGGG"""
        )
        expected_result_4 = dedent(
            """>3"""  
        )
        expected_result_5 = dedent(
            """AAATTTGGG"""
        )
        expected_result_6 = dedent(
            """>4"""  
        )
        expected_result_7 = dedent(
            """AAATTTGGG"""
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
            call(expected_result_3),
            call(expected_result_4),
            call(expected_result_5),
            call(expected_result_6),
            call(expected_result_7),
        ]

    @patch("builtins.print")
    def test_dna_threader_alias1(self, mocked_print):
        expected_result_0 = dedent(
            """>1"""
        )
        expected_result_1 = dedent(
            """AAAGGG---"""
        )
        expected_result_2 = dedent(
            """>2"""  
        )
        expected_result_3 = dedent(
            """AAATTTGGG"""
        )
        expected_result_4 = dedent(
            """>3"""  
        )
        expected_result_5 = dedent(
            """AAATTTGGG"""
        )
        expected_result_6 = dedent(
            """>4"""  
        )
        expected_result_7 = dedent(
            """AAATTTGGG"""
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
            call(expected_result_3),
            call(expected_result_4),
            call(expected_result_5),
            call(expected_result_6),
            call(expected_result_7),
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
            """>1"""
        )
        expected_result_1 = dedent(
            """AAAGGG---"""
        )
        expected_result_2 = dedent(
            """>2"""  
        )
        expected_result_3 = dedent(
            """AAATTTGGG"""
        )
        expected_result_4 = dedent(
            """>3"""  
        )
        expected_result_5 = dedent(
            """AAATTTGGG"""
        )
        expected_result_6 = dedent(
            """>4"""  
        )
        expected_result_7 = dedent(
            """AAATTTGGG"""
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
            call(expected_result_3),
            call(expected_result_4),
            call(expected_result_5),
            call(expected_result_6),
            call(expected_result_7),
        ]

    @patch("builtins.print")
    def test_dna_threader_alias_stop_codon_true(self, mocked_print):
        expected_result_0 = dedent(
            """>1"""
        )
        expected_result_1 = dedent(
            """AAAGGG---"""
        )
        expected_result_2 = dedent(
            """>2"""  
        )
        expected_result_3 = dedent(
            """AAATTTGGG"""
        )
        expected_result_4 = dedent(
            """>3"""  
        )
        expected_result_5 = dedent(
            """AAATTTGGG"""
        )
        expected_result_6 = dedent(
            """>4"""  
        )
        expected_result_7 = dedent(
            """AAATTT---"""
        )
        testargs = [
            "phykit",
            "p2n",
            "-p",
            f"{here.parent.parent.parent}/sample_files/test_alignment.prot.faa",
            "-n",
            f"{here.parent.parent.parent}/sample_files/test.nucl.fna",
            "-s",
            "0",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
            call(expected_result_2),
            call(expected_result_3),
            call(expected_result_4),
            call(expected_result_5),
            call(expected_result_6),
            call(expected_result_7),
        ]

    @patch("builtins.print")
    def test_dna_threader_alias_stop_codon_true(self, mocked_print):
        expected_result_0 = dedent(
            """>1"""
        )
        expected_result_1 = dedent(
            """AAAGGG---"""
        )
        expected_result_2 = dedent(
            """>2"""  
        )
        expected_result_3 = dedent(
            """AAATTTGGG"""
        )
        expected_result_4 = dedent(
            """>3"""  
        )
        expected_result_5 = dedent(
            """AAATTTGGG"""
        )
        expected_result_6 = dedent(
            """>4"""  
        )
        expected_result_7 = dedent(
            """AAATTTGGG"""
        )
        testargs = [
            "phykit",
            "p2n",
            "-p",
            f"{here.parent.parent.parent}/sample_files/test_alignment.prot.faa",
            "-n",
            f"{here.parent.parent.parent}/sample_files/test.nucl.fna",
            "-s",
            "1",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
            call(expected_result_2),
            call(expected_result_3),
            call(expected_result_4),
            call(expected_result_5),
            call(expected_result_6),
            call(expected_result_7),
        ]

    @patch("builtins.print")
    def test_dna_threader_ambig_char_false(self, mocked_print):
        expected_result_0 = dedent(
            """>1"""
        )
        expected_result_1 = dedent(
            """AAAGGG---"""
        )
        expected_result_2 = dedent(
            """>2"""  
        )
        expected_result_3 = dedent(
            """AAATTTGGG"""
        )
        expected_result_4 = dedent(
            """>3"""  
        )
        expected_result_5 = dedent(
            """AAATTT---"""
        )
        expected_result_6 = dedent(
            """>4"""  
        )
        expected_result_7 = dedent(
            """AAATTT---"""
        )
        testargs = [
            "phykit",
            "p2n",
            "-p",
            f"{here.parent.parent.parent}/sample_files/test_alignment.prot.ambig.faa",
            "-n",
            f"{here.parent.parent.parent}/sample_files/test.nucl.fna",
            "-s",
            "F",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
            call(expected_result_2),
            call(expected_result_3),
            call(expected_result_4),
            call(expected_result_5),
            call(expected_result_6),
            call(expected_result_7),
        ]

    @patch("builtins.print")
    def test_dna_threader_ambig_char_true(self, mocked_print):
        expected_result_0 = dedent(
            """>1"""
        )
        expected_result_1 = dedent(
            """AAAGGG---"""
        )
        expected_result_2 = dedent(
            """>2"""  
        )
        expected_result_3 = dedent(
            """AAATTTGGG"""
        )
        expected_result_4 = dedent(
            """>3"""  
        )
        expected_result_5 = dedent(
            """AAATTTGGG"""
        )
        expected_result_6 = dedent(
            """>4"""  
        )
        expected_result_7 = dedent(
            """AAATTTGGG"""
        )
        testargs = [
            "phykit",
            "p2n",
            "-p",
            f"{here.parent.parent.parent}/sample_files/test_alignment.prot.ambig.faa",
            "-n",
            f"{here.parent.parent.parent}/sample_files/test.nucl.fna",
            "-s",
            "T",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
            call(expected_result_2),
            call(expected_result_3),
            call(expected_result_4),
            call(expected_result_5),
            call(expected_result_6),
            call(expected_result_7),
        ]

    @patch("builtins.print")
    def test_dna_threader_aa_file_not_found(self, mocked_print):
        testargs = [
            "phykit",
            "p2n",
            "-p",
            "not_real",
            "-n",
            f"{here.parent.parent.parent}/sample_files/test.nucl.fna",
            "-s",
            "T",
        ]
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_dna_threader_nucl_file_not_found(self, mocked_print):
        testargs = [
            "phykit",
            "p2n",
            "-p",
            f"{here.parent.parent.parent}/sample_files/test_alignment.prot.ambig.faa",
            "-n",
            "not_real",
            "-s",
            "T",
        ]
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2