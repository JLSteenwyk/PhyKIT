import pytest
from argparse import Namespace
from math import isclose

from phykit.services.alignment.composition_per_taxon import CompositionPerTaxon


@pytest.fixture
def args():
    kwargs = dict(alignment="/some/path/to/file.fa")
    return Namespace(**kwargs)


class TestCompositionPerTaxon(object):
    def test_init_sets_alignment_file_path(self, args):
        composition = CompositionPerTaxon(args)
        assert composition.alignment_file_path == args.alignment

    def test_composition_per_taxon(self, alignment_simple, args):
        composition = CompositionPerTaxon(args)
        symbols, rows = composition.calculate_composition_per_taxon(
            alignment_simple, is_protein=False
        )

        assert symbols == ["A", "C", "G", "T"]
        by_taxon = {taxon: vals for taxon, vals in rows}

        assert isclose(by_taxon["1"][0], 0.4, rel_tol=0.001)  # A
        assert isclose(by_taxon["1"][1], 0.0, rel_tol=0.001)  # C
        assert isclose(by_taxon["1"][2], 0.2, rel_tol=0.001)  # G
        assert isclose(by_taxon["1"][3], 0.4, rel_tol=0.001)  # T

        assert isclose(by_taxon["4"][0], 0.6, rel_tol=0.001)  # A
        assert isclose(by_taxon["4"][1], 0.0, rel_tol=0.001)  # C
        assert isclose(by_taxon["4"][2], 0.2, rel_tol=0.001)  # G
        assert isclose(by_taxon["4"][3], 0.2, rel_tol=0.001)  # T
