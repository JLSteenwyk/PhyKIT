import pytest
from argparse import Namespace
from Bio import Phylo
from math import isclose

from phykit.services.tree.print_tree import PrintTree


@pytest.fixture
def args():
    kwargs = dict(tree="/some/path/to/file.tre", remove=None)
    return Namespace(**kwargs)


class TestPrintTree(object):
    def test_init_sets_tree_file_path(self, args):
        t = PrintTree(args)
        assert t.tree_file_path == args.tree
        assert t.output_file_path is None

    def test_read_file_reads_tree_file_path(self, mocker, args):
        # Mock the cached tree read method instead of Phylo.read
        mock_cached_read = mocker.patch("phykit.services.tree.base.Tree._cached_tree_read")
        mock_get_hash = mocker.patch("phykit.services.tree.base.Tree._get_file_hash", return_value="test_hash")

        t = PrintTree(args)
        t.read_tree_file()

        # Verify the cached read was called with the correct parameters
        mock_get_hash.assert_called_with(args.tree)
        mock_cached_read.assert_called_with(args.tree, "newick", "test_hash")
