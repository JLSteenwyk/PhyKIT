import pytest
from Bio import Phylo
from math import isclose

from phykit.services.tree.lb_score import LBScore

class TestLBScore(object):
    def test_init_sets_tree_file_path(self, args):
        t = LBScore(args)
        assert t.tree_file_path == args.tree
        assert t.output_file_path is None

    def test_read_file_reads_tree_file_path(self, mocker, args):
        mock_read = mocker.patch("phykit.services.tree.base.Phylo.read")
        t = LBScore(args)
        t.read_file()
        mock_read.assert_called_with(args.tree, "newick")

    def test_calculate_lb_score_zero_branch_len(self, tree_zero_branch_length, args):
        t = LBScore(args)
        res = t.calculate_lb_score(tree_zero_branch_length)
        assert res is None

    def test_calculate_treeness(self, tree_simple, args):
        t = LBScore(args)
        mean, median, twenty_fifth, seventy_fifth, standard_deviation, variance = t.calculate_lb_score(tree_simple)
        assert isinstance(mean, float)
        assert isinstance(median, float)
        assert isinstance(twenty_fifth, float)
        assert isinstance(seventy_fifth, float)
        assert isinstance(standard_deviation, float)
        assert isinstance(variance, float)
        assert isclose(mean, -12.50000000000002, rel_tol=0.001)
        assert isclose(median, -27.80598423286591, rel_tol=0.001)
        assert isclose(twenty_fifth, -31.04918307557076, rel_tol=0.001)
        assert isclose(seventy_fifth, -12.903859858133497, rel_tol=0.001)
        assert isclose(standard_deviation, 35.26687859163367, rel_tol=0.001)
        assert isclose(variance, 1243.7527255970294, rel_tol=0.001)



