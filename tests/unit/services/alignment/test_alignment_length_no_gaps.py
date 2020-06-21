import pytest
from argparse import Namespace
from math import isclose

from phykit.services.alignment.alignment_length_no_gaps import AlignmentLengthNoGaps


@pytest.fixture
def args():
    kwargs = dict(alignment="/some/path/to/file.fa")
    return Namespace(**kwargs)


class TestAlignmentLengthNoGaps(object):
    def test_init_sets_alignment_file_path(self, args):
        aln = AlignmentLengthNoGaps(args)
        assert aln.alignment_file_path == args.alignment
        assert aln.output_file_path is None

    def test_alignment_length_no_gaps(self, alignment_simple, args):
        aln = AlignmentLengthNoGaps(args)
        (
            aln_len_no_gaps,
            aln_len,
            aln_len_no_gaps_per,
        ) = aln.calculate_alignment_length_no_gaps(alignment_simple)
        assert isinstance(aln_len_no_gaps, int)
        assert isinstance(aln_len, int)
        assert isinstance(aln_len_no_gaps_per, float)
        assert aln_len_no_gaps == 3
        assert aln_len == 6
        assert isclose(aln_len_no_gaps_per, 50.0, rel_tol=0.001)
