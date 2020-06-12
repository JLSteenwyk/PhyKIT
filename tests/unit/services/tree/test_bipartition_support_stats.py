import pytest
from Bio import Phylo
from math import isclose

from phykit.services.tree.bipartition_support_stats import BipartitionSupportStats

class TestBipartitionSupportStats(object):
    def test_init_sets_tree_file_path(self, args):
        t = BipartitionSupportStats(args)
        assert t.tree_file_path == args.tree
        assert t.output_file_path is None

    def test_read_file_reads_tree_file_path(self, mocker, args):
        mock_read = mocker.patch("phykit.services.tree.base.Phylo.read")
        t = BipartitionSupportStats(args)
        t.read_tree_file()
        mock_read.assert_called_with(args.tree, "newick")

    def test_calculate_bipartition_support_stats(self, small_aspergillus_tree, args):
        t = BipartitionSupportStats(args)
        bs_vals, stats = t.calculate_bipartition_support_stats(small_aspergillus_tree)
        assert isinstance(stats['mean'], float)
        assert isinstance(stats['median'], (int or float))
        assert isinstance(stats['twenty_fifth'], float)
        assert isinstance(stats['seventy_fifth'], float)
        assert isinstance(stats['standard_deviation'], float)
        assert isinstance(stats['variance'], float)
        assert isclose(stats['mean'], 95.71428571428571, rel_tol=0.001)
        assert isclose(stats['median'], 100, rel_tol=0.001)
        assert isclose(stats['twenty_fifth'], 92.5, rel_tol=0.001)
        assert isclose(stats['seventy_fifth'], 100.0, rel_tol=0.001)
        assert isclose(stats['standard_deviation'], 7.319250547113999, rel_tol=0.001)
        assert isclose(stats['variance'], 53.57142857142857, rel_tol=0.001)
        assert isclose(stats['minimum'], 85, rel_tol=0.001)
        assert isclose(stats['maximum'], 100, rel_tol=0.001)