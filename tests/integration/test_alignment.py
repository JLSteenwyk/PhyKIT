import pytest
import sys
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestAlignment(object):
    @patch("builtins.print")
    def test_alignment_length(self, mocked_print):
        expected_result = 6
        testargs = [
            "phykit",
            "alignment_length",
            f"{here.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_alignment_length_no_gaps(self, mocked_print):
        expected_result = "3\t6\t50.0"
        testargs = [
            "phykit",
            "alignment_length_no_gaps",
            f"{here.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_create_concatenation_matrix(self, mocked_print):
        prefix = "output/create_concat_matrix"
        expected_start_message = dedent(
            """
			--------------------
			| General features |
			--------------------
			Total number of taxa: 12
			Total number of alignments: 3
			
			
			----------------
			| Output files |
			----------------
			Partition file output: output/create_concat_matrix.partition
			Concatenated fasta output: output/create_concat_matrix.fa
			Occupancy report: output/create_concat_matrix.occupancy
			"""
        )
        testargs = [
            "phykit",
            "create_concatenation_matrix",
            "-a",
            f"{here.parent.parent}/sample_files/alignment_list_for_create_concat_matrix.txt",
            "-p",
            prefix,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call(expected_start_message),
            call("Complete!\n"),
        ]

        with open(f"{here.parent}/expected/concat_matrix.fa", "r") as expected_fa, open(
            f"{here.parent}/expected/concat_matrix.occupancy", "r"
        ) as expected_occupency, open(
            f"{here.parent}/expected/concat_matrix.partition"
        ) as expected_partition:
            expected_fa_content = expected_fa.read()
            expected_occupency_content = expected_occupency.read()
            expected_partition_content = expected_partition.read()

        with open(f"{prefix}.fa", "r") as out_fa, open(
            f"{prefix}.occupancy", "r"
        ) as out_occupency, open(f"{prefix}.partition", "r") as out_partition:
            out_fa_content = out_fa.read()
            out_occupency_content = out_occupency.read()
            out_partition_content = out_partition.read()

        assert expected_fa_content == out_fa_content
        assert expected_occupency_content == out_occupency_content
        assert expected_partition_content == out_partition_content

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
            f"{here.parent.parent}/sample_files/EOG091N44MS.fa.mafft",
            "-n",
            f"{here.parent.parent}/sample_files/EOG091N44MS.fa",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
        ]


    @patch("builtins.print")
    def test_gc_content(self, mocked_print):
        expected_result = "0.22727272727272727"
        testargs = [
            "phykit",
            "gc_content",
            f"{here.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]


    @patch("builtins.print")
    def test_pairwise_identity(self, mocked_print):
        testargs = [
            "phykit",
            "pairwise_identity",
            f"{here.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call("mean: 0.4666666666666667"),
            call("median: 0.5"),
            call("25th percentile: 0.3333333333333333"),
            call("75th percentile: 0.625"),
            call("minimum: 0.16666666666666666"),
            call("maximum: 0.8333333333333334"),
            call("standard deviation: 0.2194268628681278"),
            call("variance: 0.048148148148148155")
        ]
    
    @patch("builtins.print")
    def test_parsimony_informative_sites(self, mocked_print):
        expected_result = "3\t6\t50.0"
        testargs = [
            "phykit",
            "parsimony_informative_sites",
            f"{here.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_rcv(self, mocked_print):
        expected_result = "0.36"
        testargs = [
            "phykit",
            "rcv",
            f"{here.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_rename_fasta_entries(self, mocked_print):
        testargs = [
            "phykit",
            "rename_fasta_entries",
            f"{here.parent.parent}/sample_files/simple.fa",
            "-i",
            f"{here.parent.parent}/sample_files/simple_fasta_idmap.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent}/expected/simple.fa.renamed.fa", "r") as expected_fa:
            expected_fa_content = expected_fa.read()

        with open(f"{here.parent.parent}/sample_files/simple.fa.renamed.fa", "r") as out_renamed:
            out_renamed_content = out_renamed.read()

        assert expected_fa_content == out_renamed_content

    @patch("builtins.print")
    def test_variable_sites(self, mocked_print):
        expected_result = "4\t6\t66.66666666666666"
        testargs = [
            "phykit",
            "variable_sites",
            f"{here.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]
