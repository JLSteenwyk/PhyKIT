import pytest
from argparse import Namespace
from Bio import Phylo
from math import isclose

from phykit.services.tree.lb_score import LBScore


@pytest.fixture
def args():
    kwargs = dict(tree="/some/path/to/file.tre", verbose=None)
    return Namespace(**kwargs)


class TestLBScore(object):
    def test_init_sets_tree_file_path(self, args):
        t = LBScore(args)
        assert t.tree_file_path == args.tree
        assert t.output_file_path is None

    def test_read_file_reads_tree_file_path(self, mocker, args):
        mock_read = mocker.patch("phykit.services.tree.base.Phylo.read")
        t = LBScore(args)
        t.read_tree_file()
        mock_read.assert_called_with(args.tree, "newick")

    # def test_calculate_lb_score_zero_branch_len(self, tree_zero_branch_length, args):
    #     t = LBScore(args)
    #     res = t.calculate_lb_score(tree_zero_branch_length)
    #     assert res is None

    def test_calculate_treeness(self, tree_simple, args):
        t = LBScore(args)
        tips, LBis = t.calculate_lb_score(tree_simple)
        expected_tips = [
            "raccoon",
            "bear",
            "sea_lion",
            "seal",
            "monkey",
            "cat",
            "weasel",
            "dog",
        ]
        expected_LBis = [
            -27.07902352846223,
            -39.283360704291205,
            -31.053612361805104,
            -31.04770664682598,
            65.67086344271493,
            12.796396820175548,
            -28.5329449372696,
            -21.470612084236517,
        ]
        assert tips == expected_tips
        for idx, value in enumerate(LBis):
            assert isclose(value, expected_LBis[idx], rel_tol=0.001)
