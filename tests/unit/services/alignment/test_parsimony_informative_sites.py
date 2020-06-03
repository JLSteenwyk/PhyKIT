import pytest
from argparse import Namespace
from math import isclose

from phykit.services.alignment.parsimony_informative_sites import ParsimonyInformative

class TestParsimonyInformative(object):
    def test_init_sets_alignment_file_path(self, args):
        pi = ParsimonyInformative(args)
        assert pi.alignment_file_path == args.alignment
        assert pi.output_file_path is None

    def test_parsimony_informative_sites(self, alignment_simple, args):
        pi = ParsimonyInformative(args)
        pi_sites, aln_len, pi_sites_per = pi.calculate_parsimony_informative_sites(alignment_simple)
        assert isinstance(pi_sites, int)
        assert isinstance(aln_len, int)
        assert isinstance(pi_sites_per, float)
        assert pi_sites == 3
        assert aln_len == 6
        assert isclose(pi_sites_per, 50.0, rel_tol=0.001)