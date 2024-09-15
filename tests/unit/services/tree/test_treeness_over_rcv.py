import pytest
from argparse import Namespace
from math import isclose

from phykit.services.tree.treeness_over_rcv import TreenessOverRCV


@pytest.fixture
def args():
    kwargs = dict(alignment="/some/path/to/file.fa", tree="/some/path/to/file.tre",)
    return Namespace(**kwargs)


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
