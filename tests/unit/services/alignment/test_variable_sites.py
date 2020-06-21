import pytest
from argparse import Namespace
from math import isclose

from phykit.services.alignment.variable_sites import VariableSites


@pytest.fixture
def args():
    kwargs = dict(alignment="/some/path/to/file.fa")
    return Namespace(**kwargs)


class TestVariableSites(object):
    def test_init_sets_alignment_file_path(self, args):
        vs = VariableSites(args)
        assert vs.alignment_file_path == args.alignment
        assert vs.output_file_path is None

    def test_variable_sites(self, alignment_simple, args):
        vs = VariableSites(args)
        var_sites, aln_len, var_sites_per = vs.calculate_variable_sites(
            alignment_simple
        )
        assert isinstance(var_sites, int)
        assert isinstance(aln_len, int)
        assert isinstance(var_sites_per, float)
        assert var_sites == 4
        assert aln_len == 6
        assert isclose(var_sites_per, 66.66666666666666, rel_tol=0.001)
