import pytest
from argparse import Namespace
from Bio import Phylo
from math import isclose

from phykit.services.tree.total_tree_length import TotalTreeLength


@pytest.fixture
def args():
    kwargs = dict(tree="/some/path/to/file.tre",)
    return Namespace(**kwargs)


class TestTotalTreeLength(object):
    def test_init_sets_tree_file_path(self, args):
        t = TotalTreeLength(args)
        assert t.tree_file_path == args.tree
        assert t.output_file_path is None

    def test_read_file_reads_tree_file_path(self, mocker, args):
        mock_read = mocker.patch("phykit.services.tree.base.Phylo.read")
        t = TotalTreeLength(args)
        t.read_tree_file()
        mock_read.assert_called_with(args.tree, "newick")

    def test_calculate_total_tree_length_zero_branch_len(
        self, tree_zero_branch_length, args
    ):
        t = TotalTreeLength(args)
        res = t.calculate_total_tree_length(tree_zero_branch_length)
        assert res == 0

    def test_calculate_calculate_total_tree_length(self, tree_simple, args):
        t = TotalTreeLength(args)
        res = t.calculate_total_tree_length(tree_simple)
        assert isinstance(res, (int, float))
        assert isclose(res, 277.27722, rel_tol=0.001)
