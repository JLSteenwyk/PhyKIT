import pytest
from argparse import Namespace
from math import isclose

from phykit.services.alignment.occupancy_per_taxon import OccupancyPerTaxon


@pytest.fixture
def args():
    kwargs = dict(alignment="/some/path/to/file.fa")
    return Namespace(**kwargs)


class TestOccupancyPerTaxon(object):
    def test_init_sets_alignment_file_path(self, args):
        occupancy = OccupancyPerTaxon(args)
        assert occupancy.alignment_file_path == args.alignment

    def test_occupancy_per_taxon(self, alignment_simple, args):
        occupancy = OccupancyPerTaxon(args)
        result = occupancy.calculate_occupancy_per_taxon(
            alignment_simple, is_protein=False
        )
        by_taxon = dict(result)

        assert isclose(by_taxon["1"], 5 / 6, rel_tol=0.001)
        assert isclose(by_taxon["2"], 4 / 6, rel_tol=0.001)
        assert isclose(by_taxon["3"], 4 / 6, rel_tol=0.001)
        assert isclose(by_taxon["4"], 5 / 6, rel_tol=0.001)
        assert isclose(by_taxon["5"], 4 / 6, rel_tol=0.001)
