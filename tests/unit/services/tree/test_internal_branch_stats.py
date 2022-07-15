import pytest
from argparse import Namespace
from Bio import Phylo
from math import isclose

from phykit.services.tree.internal_branch_stats import InternalBranchStats


@pytest.fixture
def args():
    kwargs = dict(tree="/some/path/to/file.tre", verbose=None)
    return Namespace(**kwargs)


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
        _, stats, _ = t.calculate_internal_branch_stats(tree_simple)
        assert isinstance(stats["mean"], float)
        assert isinstance(stats["median"], float)
        assert isinstance(stats["twenty_fifth"], float)
        assert isinstance(stats["seventy_fifth"], float)
        assert isinstance(stats["standard_deviation"], float)
        assert isinstance(stats["variance"], float)
        assert isinstance(stats["minimum"], float)
        assert isinstance(stats["maximum"], float)
        assert isclose(stats["mean"], 6.987232, rel_tol=0.001)
        assert isclose(stats["median"], 3.87382, rel_tol=0.001)
        assert isclose(stats["twenty_fifth"], 2.0946, rel_tol=0.001)
        assert isclose(stats["seventy_fifth"], 7.52973, rel_tol=0.001)
        assert isclose(stats["standard_deviation"], 8.011401268758792, rel_tol=0.001)
        assert isclose(stats["variance"], 64.18255028907, rel_tol=0.001)
        assert isclose(stats["minimum"], 0.846, rel_tol=0.001)
        assert isclose(stats["maximum"], 20.59201, rel_tol=0.001)
