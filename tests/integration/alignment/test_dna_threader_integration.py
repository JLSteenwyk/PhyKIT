import pytest
import sys
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestDNAThreader(object):
    # @patch("builtins.print")
    # def test_dna_threader0(self, mocked_print):
    #     expected_result_0 = dedent(
    #         """>200_S38"""
    #     )
    #     expected_result_1 = dedent(
    #         """atggctgacatcctcacgcagctccagacttgcctggatcagcttgcaacacaattctacgcaacacttggttatctcacaacataccacgacaatgcccccacaacaccaccacca------cccgacgcagcaccagccctagcaaagatcaccaagaactcatcatcaccgccagtcccagcagccatcgcaaataaagtggggggtgcagctgctgttgcgggcaatgcatcaccc---caggcgcctcctcaacaa---------gcg---------------------------------------------ccagactcgcccagc---------agccgg---cgggagcttgcgcgcgatctcattatcaaagaacagcagatcgagtaccttatctccgtgcttcccgggattggcgcctctgaggctgaacaagaaaccagaatccaggacctggagaccgagcttagagacgtcgagaaggagcgcgctgcgaaagtgcgggagttgaaaaagttgaggactcggttggaggatgttcttggcgctgtcgctgtgggtatccacggggatggttactctcaaaactga"""
    #     )
    #     expected_result_2 = dedent(
    #         """>203_S40"""
    #     )
    #     expected_result_3 = dedent(
    #         """atggctgacatcctcacgcagctccagacttgcctggatcagcttgcaacacaattctacgcaacacttggttatctcacaacataccacgacaatgcccccacaacaccaccacca------cccgacgcagcaccagccctagcaaagatcaccaagaactcatcatcaccaccagtcccagcagccatcgcaaataaagtggggggtgcagctgctgttgcgggcaatgcatcaccc---caggcgcctcctcaacaa---------gcg---------------------------------------------ccagactcgcccagc---------agccgg---cgggagcttgcgcgcgatctcattatcaaagaacagcagatcgagtaccttatctccgtgcttcccgggattggcgcctctgaggctgaacaagaaaccagaatccaggacctggagaccgagcttagagacgtcgagaaggagcgcgctgcgaaagtgcgggagttgaaaaagttgaggactcggttggaggatgttcttggcgctgtcgctgtgggtatccacggggatggttactctcaaaactga"""
    #     )

    #     testargs = [
    #         "phykit",
    #         "thread_dna",
    #         "-p",
    #         f"{here.parent.parent.parent}/sample_files/EOG091N44MS.fa.mafft",
    #         "-n",
    #         f"{here.parent.parent.parent}/sample_files/EOG091N44MS.fa",
    #     ]

    #     with patch.object(sys, "argv", testargs):
    #         Phykit()
    #     assert mocked_print.mock_calls == [
    #         call(expected_result_0),
    #         call(expected_result_1),
    #         call(expected_result_2),
    #         call(expected_result_3),
    #     ]

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
    def test_dna_threader_alias_stop_codon_true1(self, mocked_print):
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
            "1"
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
    def test_dna_threader_alias_stop_codon_false(self, mocked_print):
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
    def test_dna_threader_alias_stop_codon_true0(self, mocked_print):
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
            """AAATTT---"""
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

    @patch("builtins.print")
    def test_dna_threader_trimmed_alignment_short_log(self, mocked_print):
        expected_result_0 = dedent(
            """>1"""
        )
        expected_result_1 = dedent(
            """AAATTT"""
        )
        expected_result_2 = dedent(
            """>2"""
        )
        expected_result_3 = dedent(
            """AAAGGG"""
        )
        expected_result_4 = dedent(
            """>3"""
        )
        expected_result_5 = dedent(
            """AAA---"""
        )
        expected_result_6 = dedent(
            """>4"""
        )
        expected_result_7 = dedent(
            """AAAGGG"""
        )
        testargs = [
            "phykit",
            "p2n",
            "-p",
            f"{here.parent.parent.parent}/sample_files/test_alignment.prot.faa.clipkit",
            "-n",
            f"{here.parent.parent.parent}/sample_files/test.nucl.fna",
            "-c",
           f"{here.parent.parent.parent}/sample_files/test_alignment.prot.faa.clipkit.log",
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
    def test_dna_threader_trimmed_alignment_long_log(self, mocked_print):
        expected_result_0 = dedent(
            """>1"""
        )
        expected_result_1 = dedent(
            """AAATTT"""
        )
        expected_result_2 = dedent(
            """>2"""  
        )
        expected_result_3 = dedent(
            """AAAGGG"""
        )
        expected_result_4 = dedent(
            """>3"""  
        )
        expected_result_5 = dedent(
            """AAA---"""
        )
        expected_result_6 = dedent(
            """>4"""  
        )
        expected_result_7 = dedent(
            """AAAGGG"""
        )
        testargs = [
            "phykit",
            "p2n",
            "-p",
            f"{here.parent.parent.parent}/sample_files/test_alignment.prot.faa.clipkit",
            "-n",
            f"{here.parent.parent.parent}/sample_files/test.nucl.fna",
            "--clipkit_log_file",
           f"{here.parent.parent.parent}/sample_files/test_alignment.prot.faa.clipkit.log",
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
    def test_dna_threader_trimmed_longer_alignment_log(self, mocked_print):
        expected_result_0 = dedent(
            """>1"""
        )
        expected_result_1 = dedent(
            """AAATTTCCCAAA------GGG---"""
            # AAAGGGTTTGGGAAACCCAAAGGG---TTT---GGG---AAACCCAAATTTGGGAAAGGGTTTGGGAAAGGG
            # 111222333444555666777888999
            # AAA===TTT======CCCAAA===---===---GGGAAACCCAAATTTGGGAAAGGGTTTGGGAAAGGG
            # MMM###LLL######GGGGGG###---###---GGGXXX

            # 1 keep Const 0.0
            # 2 trim nConst,nPI 0.0
            # 3 keep nConst,nPI 0.25
            # 4 trim Const 0.0
            # 5 trim Const 0.0
            # 6 keep Const 0.0
            # 7 keep Const 0.0
            # 8 trim Const 0.0
            # 9 keep Const 0.0
            # 10 trim Const 0.0
            # 11 keep Const 0.0
            # 12 keep Const 0.0
            # 13 keep Const 0.0
        )
        expected_result_2 = dedent(
            """>2"""  
        )
        expected_result_3 = dedent(
            """AAAGGGTTTGGGAAAGGGAAAGGG"""
        )
        expected_result_4 = dedent(
            """>3"""  
        )
        expected_result_5 = dedent(
            """AAAGGGTTTGGGAAAGGGAAAGGG"""
        )
        expected_result_6 = dedent(
            """>4"""  
        )
        expected_result_7 = dedent(
            """AAAGGGTTTGGGAAAGGGAAAGGG"""
        )
        testargs = [
            "phykit",
            "p2n",
            "-p",
            f"{here.parent.parent.parent}/sample_files/test_alignment_longer.prot.faa.clipkit",
            "-n",
            f"{here.parent.parent.parent}/sample_files/test_longer.nucl.fna",
            "--clipkit_log_file",
           f"{here.parent.parent.parent}/sample_files/test_alignment_longer.prot.faa.clipkit.log",
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
