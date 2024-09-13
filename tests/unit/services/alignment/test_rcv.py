import pytest
from argparse import Namespace
from math import isclose

from phykit.services.alignment.rcv import RelativeCompositionVariability


@pytest.fixture
def args():
    kwargs = dict(alignment="/some/path/to/file.fa")
    return Namespace(**kwargs)


class TestRelativeCompositionVariability(object):
    def test_init_sets_alignment_file_path(self, args):
        rcv = RelativeCompositionVariability(args)
        assert rcv.alignment_file_path == args.alignment
        assert rcv.output_file_path is None

    def test_relative_composition_variability(self, mocker, alignment_simple, args):
        mocker.patch(
            "phykit.services.alignment.rcv.RelativeCompositionVariability.get_alignment_and_format",
            return_value=(alignment_simple, "fa", True),
        )
        rcv = RelativeCompositionVariability(args)
        relative_composition_variability = rcv.calculate_rcv()
        assert isinstance(relative_composition_variability, float)
        assert isclose(relative_composition_variability, 0.36, rel_tol=0.001)
