import pytest
from Bio import Phylo
from math import isclose

from phykit.services.tree.patristic_distances import PatristicDistances

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
        mean, median, twenty_fifth, seventy_fifth, minimum, maximum, standard_deviation, variance, patristic_distances, combos = t.calculate_patristic_distances(tree_simple)
        assert isinstance(mean, float)
        assert isinstance(median, float)
        assert isinstance(twenty_fifth, float)
        assert isinstance(seventy_fifth, float)
        assert isinstance(standard_deviation, float)
        assert isinstance(variance, float)
        assert isinstance(minimum, float)
        assert isinstance(maximum, float)
        assert isclose(mean, 76.19737857142857, rel_tol=0.001)
        assert isclose(median, 49.588789999999996, rel_tol=0.001)
        assert isclose(twenty_fifth, 40.50536, rel_tol=0.001)
        assert isclose(seventy_fifth, 108.13853, rel_tol=0.001)
        assert isclose(standard_deviation, 45.46979239234539, rel_tol=0.001)
        assert isclose(variance, 2067.5020202029905, rel_tol=0.001)
        assert isclose(minimum, 24.0, rel_tol=0.001)
        assert isclose(maximum, 152.88127, rel_tol=0.001)