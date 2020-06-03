import pytest
from Bio import Phylo
from math import isclose

from phykit.services.tree.treeness_over_rcv import TreenessOverRCV


class TestTreenessOverRCV(object):
    def test_init_sets_tree_file_path(self, args):
        t = TreenessOverRCV(args)
        assert t.tree_file_path == args.tree
        assert t.alignment_file_path == args.alignment        
        assert t.output_file_path is None

    def test_read_file_reads_tree_file_path(self, mocker, args):
        mock_read = mocker.patch("phykit.services.tree.base.Phylo.read")
        t = TreenessOverRCV(args)
        t.read_tree_file()
        mock_read.assert_called_with(args.tree, "newick")

    def test_calculate_treeness_over_rcv(self, args):
        t = TreenessOverRCV(args)
        treeness = 0.12599722400563595
        rcv = 0.36
        res = t.calculate_treeness_over_rcv(treeness, rcv)
        assert isinstance(res, float)
        assert isclose(res, 0.3499922889045443, rel_tol=0.001)