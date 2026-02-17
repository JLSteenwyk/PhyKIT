import pytest
from argparse import Namespace

from phykit.services.alignment.mask_alignment import MaskAlignment


@pytest.fixture
def args():
    kwargs = dict(
        alignment="/some/path/to/file.fa",
        max_gap=1.0,
        min_occupancy=0.0,
        max_entropy=None,
    )
    return Namespace(**kwargs)


class TestMaskAlignment(object):
    def test_init_sets_alignment_file_path(self, args):
        masker = MaskAlignment(args)
        assert masker.alignment_file_path == args.alignment

    def test_keep_mask_by_gap_and_occupancy(self, alignment_simple, args):
        args.max_gap = 0.3
        args.min_occupancy = 0.8
        masker = MaskAlignment(args)
        keep_mask = masker.calculate_keep_mask(alignment_simple, is_protein=False)
        assert keep_mask.tolist() == [True, False, True, False, True, True]

    def test_keep_mask_by_entropy(self, alignment_simple, args):
        args.max_entropy = 0.5
        masker = MaskAlignment(args)
        keep_mask = masker.calculate_keep_mask(alignment_simple, is_protein=False)
        assert keep_mask.tolist() == [True, False, False, True, False, False]
