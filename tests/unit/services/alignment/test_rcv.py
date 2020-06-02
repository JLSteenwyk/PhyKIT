import pytest
from argparse import Namespace
from math import isclose

from phykit.services.alignment.rcv import RelativeCompositionVariability

class TestRelativeCompositionVariability(object):
    def test_init_sets_alignment_file_path(self, args):
        rcv = RelativeCompositionVariability(args)
        assert rcv.alignment_file_path == args.alignment
        assert rcv.output_file_path is None

    def test_relative_composition_variability(self, alignment_simple, args):
        rcv = RelativeCompositionVariability(args)
        aln_len = 6
        relative_composition_variability = rcv.calculate_rcv(alignment_simple, aln_len)
        assert isinstance(relative_composition_variability, float)
        assert isclose(relative_composition_variability, 0.36, rel_tol=0.001)