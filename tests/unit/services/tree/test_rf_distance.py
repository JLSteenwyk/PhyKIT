import pytest
from argparse import Namespace
from Bio import Phylo
from math import isclose

from phykit.services.tree.rf_distance import RobinsonFouldsDistance


@pytest.fixture
def args():
    kwargs = dict(
        tree_zero="/some/path/to/file.tre", tree_one="/some/path/to/file.tre",
    )
    return Namespace(**kwargs)


class TestRobinsonFouldsDistance(object):
    def test_init_sets_tree_file_path(self, args):
        rf = RobinsonFouldsDistance(args)
        assert rf.tree_file_path == args.tree_zero
        assert rf.tree1_file_path == args.tree_one
        assert rf.output_file_path is None

    def test_read_file_reads_tree_file_path(self, mocker, args):
        # Mock the cached tree read method instead of Phylo.read
        mock_cached_read = mocker.patch("phykit.services.tree.base.Tree._cached_tree_read")
        mock_get_hash = mocker.patch("phykit.services.tree.base.Tree._get_file_hash", return_value="test_hash")

        rf = RobinsonFouldsDistance(args)
        rf.read_tree_file()

        # Verify the cached read was called with the correct parameters
        mock_get_hash.assert_called_with(args.tree_zero)
        mock_cached_read.assert_called_with(args.tree_zero, "newick", "test_hash")

    def test_calculate_patristic_distances(self, tree_simple, tree_simple_other, args):
        rf = RobinsonFouldsDistance(args)
        plain_rf, normalized_rf = rf.calculate_robinson_foulds_distance(
            tree_simple, tree_simple_other
        )
        assert isinstance(plain_rf, int)
        assert isinstance(normalized_rf, float)
        assert plain_rf == 8
        assert normalized_rf == 0.8
