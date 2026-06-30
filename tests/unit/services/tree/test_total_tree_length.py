import pytest
import subprocess
import sys
from argparse import Namespace
from math import isclose

from phykit.services.tree.total_tree_length import TotalTreeLength


def test_module_import_does_not_import_json_or_typing():
    code = """
import sys
import phykit.services.tree.total_tree_length as module
assert callable(module.print_json)
assert "json" not in sys.modules
assert "typing" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


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
        # Mock the cached tree read method instead of Phylo.read
        mock_cached_read = mocker.patch("phykit.services.tree.base.Tree._cached_tree_read")
        mock_get_hash = mocker.patch("phykit.services.tree.base.Tree._get_file_hash", return_value="test_hash")

        t = TotalTreeLength(args)
        t.read_tree_file()

        # Verify the cached read was called with the correct parameters
        mock_get_hash.assert_called_with(args.tree)
        mock_cached_read.assert_called_with(args.tree, "newick", "test_hash")

    def test_run_uses_unmodified_tree_read(self, mocker, args):
        tree = object()
        t = TotalTreeLength(args)
        read_tree = mocker.patch.object(
            t,
            "read_tree_file_unmodified",
            return_value=tree,
        )
        mocker.patch.object(t, "calculate_total_tree_length", return_value=12.34567)
        mocked_print = mocker.patch("builtins.print")

        t.run()

        read_tree.assert_called_once_with()
        mocked_print.assert_called_once_with(12.3457)

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
