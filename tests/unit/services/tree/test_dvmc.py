import pytest
from argparse import Namespace
from Bio import Phylo
from math import isclose

from phykit.services.tree.dvmc import DVMC


@pytest.fixture
def args():
    kwargs = dict(tree="/some/path/to/file.tre")
    return Namespace(**kwargs)


class TestDVMC(object):
    def test_init_sets_tree_file_path(self, args):
        d = DVMC(args)
        assert d.output_file_path is None

    def test_read_file_reads_tree_file_path(self, mocker, args):
        mock_read = mocker.patch("phykit.services.tree.base.Phylo.read")
        d = DVMC(args)
        d.read_tree_file()
        mock_read.assert_called_with(args.tree, "newick")

