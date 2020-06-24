import pytest
from argparse import Namespace

from phykit.services.alignment.create_concatenation_matrix import CreateConcatenationMatrix

@pytest.fixture
def args():
    kwargs = dict(
        alignment_list="/some/path/to/file",
        prefix="some_str")
    return Namespace(**kwargs)

class TestCreateConcatenationMatrix(object):
    def test_init_sets_alignment_list_path(self, args):
        ccm = CreateConcatenationMatrix(args)
        assert ccm.alignment_list == args.alignment_list
        assert ccm.prefix == args.prefix
        assert ccm.output_file_path is None
    
    def test_taxa_names_acquisition(self, alignments, args):
        expected_taxa = [
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

        ccm = CreateConcatenationMatrix(args)
        taxa = ccm.get_taxa_names(alignments)

        assert isinstance(taxa, list)
        assert taxa == expected_taxa