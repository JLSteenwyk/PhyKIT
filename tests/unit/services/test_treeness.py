import pytest
from argparse import Namespace
from Bio import Phylo
from math import isclose
from pathlib import Path

from phykit.services.treeness import Treeness

here = Path(__file__)


@pytest.fixture
def args():
    kwargs = dict(tree="/some/path/to/file.tre")
    return Namespace(**kwargs)


@pytest.fixture
def tree_zero_branch_length(mocker):
    return Phylo.read(
        f"{here.parent.parent.parent}/sample_files/tree_zero_branch_length.tre",
        "newick",
    )


@pytest.fixture
def tree_simple(mocker):
    return Phylo.read(
        f"{here.parent.parent.parent}/sample_files/tree_simple.tre", "newick",
    )


class TestTreeness(object):
    def test_process_args_sets_tree_file_path(self, args):
        t = Treeness(args)
        assert t.tree_file_path == args.tree

    def test_read_file_reads_tree_file_path(self, mocker, args):
        mock_read = mocker.patch("phykit.services.treeness.Phylo.read")
        t = Treeness(args)
        t.read_file()
        mock_read.assert_called_with(args.tree, "newick")

    def test_calculate_treeness_zero_branch_len(self, tree_zero_branch_length, args):
        t = Treeness(args)
        res = t.calculate_treeness(tree_zero_branch_length)
        assert res is None

    def test_calculate_treeness(self, tree_simple, args):
        t = Treeness(args)
        res = t.calculate_treeness(tree_simple)
        assert isinstance(res, float)
        assert isclose(res, 0.12599722400563595, rel_tol=0.001)
