import pytest
from Bio import Phylo
from math import isclose

from phykit.services.tree.dvmc import DVMC


class TestDVMC(object):
    def test_init_sets_tree_file_path(self, args):
        d = DVMC(args)
        assert d.outgroup_taxa_file_path == args.root
        assert d.output_file_path is None

    def test_read_file_reads_tree_file_path(self, mocker, args):
        mock_read = mocker.patch("phykit.services.tree.base.Phylo.read")
        d = DVMC(args)
        d.read_tree_file()
        mock_read.assert_called_with(args.tree, "newick")

    def test_calculate_dvmc_zero_branch_len(self, tree_zero_branch_length, args):
        d = DVMC(args)
        res = d.calculate_dvmc(tree_zero_branch_length, args.root)
        assert res == 0.0

    def test_calculate_dvmc(self, tree_simple, tree_simple_outgroup, args):
        d = DVMC(args)
        res = d.calculate_dvmc(tree_simple, tree_simple_outgroup)
        assert isinstance(res, float)
        assert isclose(res, 42.80162365633575, rel_tol=0.001)