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
        expected_result = "6"
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
            "output/create_concat_matrix",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_start_message),
            call("Complete!\n"),
        ]
