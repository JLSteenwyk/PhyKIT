import pytest
import sys
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestCreateConcatenationMatrix(object):
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
            f"{here.parent.parent.parent}/sample_files/alignment_list_for_create_concat_matrix.txt",
            "-p",
            prefix,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call(expected_start_message),
            call("Complete!\n"),
        ]

        with open(f"{here.parent.parent}/expected/concat_matrix.fa", "r") as expected_fa, open(
            f"{here.parent.parent}/expected/concat_matrix.occupancy", "r"
        ) as expected_occupency, open(
            f"{here.parent.parent}/expected/concat_matrix.partition"
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