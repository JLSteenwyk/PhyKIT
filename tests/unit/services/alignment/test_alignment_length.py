import pytest
from argparse import Namespace
from math import isclose

from phykit.services.alignment.alignment_length import AlignmentLength
from phykit.services.alignment.base import Alignment


@pytest.fixture
def args():
    kwargs = dict(alignment="/some/path/to/file.fa")
    return Namespace(**kwargs)


class TestAlignmentLength(object):
    def test_init_sets_alignment_file_path(self, args):
        aln = AlignmentLength(args)
        assert aln.alignment_file_path == args.alignment
        assert aln.output_file_path is None

    def test_alignment_length(self, alignment_simple, args):
        aln = AlignmentLength(args)
        res = aln.calculate_alignment_length(alignment_simple)
        assert isinstance(res, int)
        assert res == 6
