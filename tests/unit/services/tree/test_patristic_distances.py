import pytest
from argparse import Namespace
from Bio import Phylo
from math import isclose

from phykit.services.tree.patristic_distances import PatristicDistances


@pytest.fixture
def args():
    kwargs = dict(tree="/some/path/to/file.tre", verbose=None)
    return Namespace(**kwargs)


class TestPatristicDistances(object):
    def test_init_sets_tree_file_path(self, args):
        t = PatristicDistances(args)
        assert t.tree_file_path == args.tree
        assert t.output_file_path is None

    def test_read_file_reads_tree_file_path(self, mocker, args):
        mock_read = mocker.patch("phykit.services.tree.base.Phylo.read")
        t = PatristicDistances(args)
        t.read_tree_file()
        mock_read.assert_called_with(args.tree, "newick")

    def test_calculate_patristic_distances(self, tree_simple, args):
        t = PatristicDistances(args)
        patristic_distances, combos, stats = t.calculate_patristic_distances(
            tree_simple
        )
        assert isinstance(stats["mean"], float)
        assert isinstance(stats["median"], float)
        assert isinstance(stats["twenty_fifth"], float)
        assert isinstance(stats["seventy_fifth"], float)
        assert isinstance(stats["standard_deviation"], float)
        assert isinstance(stats["variance"], float)
        assert isinstance(stats["minimum"], float)
        assert isinstance(stats["maximum"], float)
        assert isclose(stats["mean"], 76.19737857142857, rel_tol=0.001)
        assert isclose(stats["median"], 49.588789999999996, rel_tol=0.001)
        assert isclose(stats["twenty_fifth"], 40.50536, rel_tol=0.001)
        assert isclose(stats["seventy_fifth"], 108.13853, rel_tol=0.001)
        assert isclose(stats["standard_deviation"], 45.46979239234539, rel_tol=0.001)
        assert isclose(stats["variance"], 2067.5020202029905, rel_tol=0.001)
        assert isclose(stats["minimum"], 24.0, rel_tol=0.001)
        assert isclose(stats["maximum"], 152.88127, rel_tol=0.001)
