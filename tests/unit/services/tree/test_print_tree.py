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
        mock_read = mocker.patch("phykit.services.tree.base.Phylo.read")
        t = PrintTree(args)
        t.read_tree_file()
        mock_read.assert_called_with(args.tree, "newick")
