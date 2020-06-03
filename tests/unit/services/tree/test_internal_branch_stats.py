import pytest
from Bio import Phylo
from math import isclose

from phykit.services.tree.internal_branch_stats import InternalBranchStats

class TestLBScore(object):
    def test_init_sets_tree_file_path(self, args):
        t = InternalBranchStats(args)
        assert t.tree_file_path == args.tree
        assert t.output_file_path is None

    def test_read_file_reads_tree_file_path(self, mocker, args):
        mock_read = mocker.patch("phykit.services.tree.base.Phylo.read")
        t = InternalBranchStats(args)
        t.read_tree_file()
        mock_read.assert_called_with(args.tree, "newick")

    def test_calculate_internal_branch_stats(self, tree_simple, args):
        t = InternalBranchStats(args)
        mean, median, twenty_fifth, seventy_fifth, minimum, maximum, standard_deviation, variance, internal_branch_lengths = t.calculate_internal_branch_stats(tree_simple)
        assert isinstance(mean, float)
        assert isinstance(median, float)
        assert isinstance(twenty_fifth, float)
        assert isinstance(seventy_fifth, float)
        assert isinstance(standard_deviation, float)
        assert isinstance(variance, float)
        assert isinstance(minimum, float)
        assert isinstance(maximum, float)
        assert isclose(mean, 6.987232, rel_tol=0.001)
        assert isclose(median, 3.87382, rel_tol=0.001)
        assert isclose(twenty_fifth, 2.0946, rel_tol=0.001)
        assert isclose(seventy_fifth, 7.52973, rel_tol=0.001)
        assert isclose(standard_deviation, 8.011401268758792, rel_tol=0.001)
        assert isclose(variance, 64.18255028907, rel_tol=0.001)
        assert isclose(minimum, 0.846, rel_tol=0.001)
        assert isclose(maximum, 20.59201, rel_tol=0.001)



