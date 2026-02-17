import pytest
from argparse import Namespace
from math import isclose

from phykit.services.alignment.alignment_entropy import AlignmentEntropy


@pytest.fixture
def args():
    kwargs = dict(alignment="/some/path/to/file.fa", verbose=False)
    return Namespace(**kwargs)


class TestAlignmentEntropy(object):
    def test_init_sets_alignment_file_path(self, args):
        entropy = AlignmentEntropy(args)
        assert entropy.alignment_file_path == args.alignment
        assert entropy.verbose is False

    def test_site_entropies(self, alignment_simple, args):
        entropy = AlignmentEntropy(args)
        entropies = entropy.calculate_site_entropies(alignment_simple, is_protein=False)

        assert len(entropies) == 6
        assert isclose(entropies[0], 0.0, rel_tol=0.001)
        assert isclose(entropies[1], 1.0, rel_tol=0.001)
        assert isclose(entropies[2], 0.97095, rel_tol=0.001)
        assert isclose(entropies[3], 0.0, rel_tol=0.001)
        assert isclose(entropies[4], 0.97095, rel_tol=0.001)
        assert isclose(entropies[5], 1.0, rel_tol=0.001)
