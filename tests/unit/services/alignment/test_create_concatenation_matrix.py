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
        assert ccm.alignment_list_path == args.alignment_list
        assert ccm.prefix == args.prefix
        assert ccm.output_file_path is None
    
    def test_taxa_names_acquisition(self, alignments, args, taxa_list):

        ccm = CreateConcatenationMatrix(args)
        taxa = ccm.get_taxa_names(alignments)

        assert isinstance(taxa, list)
        assert sorted(taxa) == sorted(taxa_list)

    def test_print_start_message(self, alignments, args, taxa_list, capfd):

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
        ccm.print_start_message(
            taxa_list,
            alignments,
            file_partition,
            fasta_output,
            file_occupancy
        )

        out, _ = capfd.readouterr()
        assert out == expected_start_message

    # def test_create_concatenation_matrix(self, args, mocker):
    #     prefix = "my_cool_prefix"
    #     alignment_list_path = "my/swaggy/path.txt"
    #     file_partition_filename = f"{prefix}.partition"
    #     fasta_output_filename = f"{prefix}.fa"
    #     file_occupancy_filename = f"{prefix}.occupancy"
    #     alignment_paths = ['path/to/file_A.fa', 'path/to/file_B.fa', 'path/to/file_C.fa']
    #     taxa_names = ['A', 'B', 'C', 'D']

    #     mocked_read_alignment_paths = mocker.path(
    #         "phykit.services.alignment.create_concatenation_matrix.CreateConcatenationMatrix.read_alignment_paths",
    #         return_value=alignment_paths
    #     )
    #     mocked_get_taxa_names = mocker.path(
    #         "phykit.services.alignment.create_concatenation_matrix.CreateConcatenationMatrix.get_taxa_names",
    #         return_value=taxa_names
    #     )
    #     mocked_print_start_message = mocker.patch(
    #         "phykit.services.alignment.create_concatenation_matrix.CreateConcatenationMatrix.print_start_message"
    #     )
    #     mocked_create_output_files = mocker.patch(
    #         "phykit.services.alignment.create_concatenation_matrix.CreateConcatenationMatrix.create_output_files"
    #     )
    #     mocked_get_list_of_taxa_and_records = mocker.patch(
    #         "phykit.services.alignment.create_concatenation_matrix.CreateConcatenationMatrix.get_list_of_taxa_and_records"
    #     )
    #     ccm = CreateConcatenationMatrix(args)
    #     res = ccm.create_concatenation_matrix(alignment_list_path, prefix)
    #     assert res is None
    #     assert mocked_print_start_message.called_with(
    #         taxa_names,
    #         alignment_paths,
    #         file_partition_filename,
    #         fasta_output_filename,
    #         file_occupancy_filename
    #     )
    #     assert mocked_get_taxa_names.called_with(
    #         alignment_paths
    #     )
    #     assert mocked_create_output_files.called_with(
    #         file_partition_filename,
    #         file_occupancy_filename
    #     )


