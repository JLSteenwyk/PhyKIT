import pytest
import sys

from argparse import Namespace
from Bio import Phylo
from math import isclose
from mock import patch, call

from phykit.phykit import Phykit
from phykit.services.tree.bipartition_support_stats import BipartitionSupportStats

@pytest.fixture
def args():
    kwargs = dict(tree="/some/path/to/file.tre", verbose=None)
    return Namespace(**kwargs)


class TestBipartitionSupportStats(object):
    def test_init_sets_tree_file_path(self, args):
        t = BipartitionSupportStats(args)
        assert t.tree_file_path == args.tree
        assert t.output_file_path is None

    def test_read_file_reads_tree_file_path(self, mocker, args):
        # Mock the cached tree read method instead of Phylo.read
        mock_cached_read = mocker.patch("phykit.services.tree.base.Tree._cached_tree_read")
        mock_get_hash = mocker.patch("phykit.services.tree.base.Tree._get_file_hash", return_value="test_hash")

        t = BipartitionSupportStats(args)
        t.read_tree_file()

        # Verify the cached read was called with the correct parameters
        mock_get_hash.assert_called_with(args.tree)
        mock_cached_read.assert_called_with(args.tree, "newick", "test_hash")

    @patch("builtins.print")
    def test_bad_file_path(self, mocked_print):
        testargs = [
            "phykit",
            "bss",
            "some/file/path",
        ]

        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type == SystemExit
        mocked_print.assert_has_calls([
            call("Please check filename and pathing"),
        ])


    # def test_calculate_bipartition_support_stats(self, small_aspergillus_tree, args):
    #     t = BipartitionSupportStats(args)
    #     bs_vals = t.get_bipartition_support_vals(small_aspergillus_tree)
    #     assert isinstance(stats["mean"], float)
    #     assert isinstance(stats["median"], (int or float))
    #     assert isinstance(stats["twenty_fifth"], float)
    #     assert isinstance(stats["seventy_fifth"], float)
    #     assert isinstance(stats["standard_deviation"], float)
    #     assert isinstance(stats["variance"], float)
    #     assert isclose(stats["mean"], 95.71428571428571, rel_tol=0.001)
    #     assert isclose(stats["median"], 100, rel_tol=0.001)
    #     assert isclose(stats["twenty_fifth"], 92.5, rel_tol=0.001)
    #     assert isclose(stats["seventy_fifth"], 100.0, rel_tol=0.001)
    #     assert isclose(stats["standard_deviation"], 7.319250547113999, rel_tol=0.001)
    #     assert isclose(stats["variance"], 53.57142857142857, rel_tol=0.001)
    #     assert isclose(stats["minimum"], 85, rel_tol=0.001)
    #     assert isclose(stats["maximum"], 100, rel_tol=0.001)
