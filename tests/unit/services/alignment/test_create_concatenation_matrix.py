from argparse import Namespace
import io
import pytest
import sys
from textwrap import dedent

from phykit.services.alignment.create_concatenation_matrix import CreateConcatenationMatrix

@pytest.fixture
def args():
    kwargs = dict(
        alignment_list="/some/path/to/file",
        prefix="some_str")
    return Namespace(**kwargs)

@pytest.fixture
def taxa_list():
    return [
        'Kpol',
        'Kpha',
        'Kbla',
        'Sdai',
        'Scas',
        'Snag',
        'Kafr',
        'Cgla',
        'Suva',
        'Skud',
        'Smik',
        'Scer'
    ]

class TestCreateConcatenationMatrix(object):
    def test_init_sets_alignment_list_path(self, args):
        ccm = CreateConcatenationMatrix(args)
        assert ccm.alignment_list == args.alignment_list
        assert ccm.prefix == args.prefix
        assert ccm.output_file_path is None
    
    def test_taxa_names_acquisition(self, alignments, args, taxa_list):

        ccm = CreateConcatenationMatrix(args)
        taxa = ccm.get_taxa_names(alignments)

        assert isinstance(taxa, list)
        assert taxa == taxa_list

    def test_start_message_print(self, alignments, args, taxa_list, capfd):

        expected_start_message = dedent(f"""
            --------------------
            | General features |
            --------------------
            Total number of taxa: 12
            Total number of alignments: 3


            ----------------
            | Output files |
            ----------------
            Partition file output: test.partition
            Concatenated fasta output: test.fa
            Occupancy report: test.occupancy

        """)

        file_partition = 'test.partition'
        fasta_output = 'test.fa'
        file_occupancy = 'test.occupancy'

        ccm = CreateConcatenationMatrix(args)
        ccm.start_message_print(
            taxa_list,
            alignments,
            file_partition,
            fasta_output,
            file_occupancy
        )

        out, _ = capfd.readouterr()
        assert out == expected_start_message

