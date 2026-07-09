import pytest
import subprocess
import sys
from argparse import Namespace
from math import isclose

from phykit.services.tree.treeness import Treeness
import phykit.services.tree.treeness as treeness_module


def test_module_import_does_not_import_json_or_typing():
    code = """
import sys
import phykit.services.tree.treeness as module
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


class TestTreeness(object):
    def test_init_sets_tree_file_path(self, args):
        t = Treeness(args)
        assert t.tree_file_path == args.tree
        assert t.output_file_path is None

    def test_read_file_reads_tree_file_path(self, mocker, args):
        # Mock the cached tree read method instead of Phylo.read
        mock_cached_read = mocker.patch("phykit.services.tree.base.Tree._cached_tree_read")
        mock_get_hash = mocker.patch("phykit.services.tree.base.Tree._get_file_hash", return_value="test_hash")

        t = Treeness(args)
        t.read_tree_file()

        # Verify the cached read was called with the correct parameters
        mock_get_hash.assert_called_with(args.tree)
        mock_cached_read.assert_called_with(args.tree, "newick", "test_hash")

    def test_calculate_treeness_zero_branch_len(self, tree_zero_branch_length, args):
        t = Treeness(args)
        res = t.calculate_treeness(tree_zero_branch_length)
        assert res is None

    def test_calculate_treeness(self, tree_simple, args):
        t = Treeness(args)
        res = t.calculate_treeness(tree_simple)
        assert isinstance(res, float)
        assert isclose(res, 0.12599722400563595, rel_tol=0.001)

    def test_process_args_defaults_json_false(self):
        parsed = Treeness(Namespace(tree="x.tre")).process_args(Namespace(tree="x.tre"))
        assert parsed["json_output"] is False

    def test_run_prints_value(self, mocker, capsys):
        t = Treeness(Namespace(tree="x.tre", json=False))
        mocker.patch.object(t, "_get_simple_newick_summary", return_value=None)
        mocker.patch.object(Treeness, "read_tree_file_unmodified", return_value=object())
        mocker.patch.object(Treeness, "calculate_treeness", return_value=0.987654)
        t.run()
        out, _ = capsys.readouterr()
        assert out.strip() == "0.9877"

    def test_run_uses_unmodified_tree_read(self, mocker):
        tree = object()
        t = Treeness(Namespace(tree="x.tre", json=False))
        mocker.patch.object(t, "_get_simple_newick_summary", return_value=None)
        read_tree = mocker.patch.object(
            t,
            "read_tree_file_unmodified",
            return_value=tree,
        )
        mocker.patch.object(t, "calculate_treeness", return_value=0.987654)
        mocker.patch("builtins.print")

        t.run()

        read_tree.assert_called_once_with()

    def test_run_json(self, mocker):
        t = Treeness(Namespace(tree="x.tre", json=True))
        mocker.patch.object(t, "_get_simple_newick_summary", return_value=None)
        mocker.patch.object(Treeness, "read_tree_file_unmodified", return_value=object())
        mocker.patch.object(Treeness, "calculate_treeness", return_value=0.11111)
        mocked_json = mocker.patch.object(treeness_module, "print_json")
        t.run()
        mocked_json.assert_called_once_with({"treeness": 0.1111})

    def test_run_uses_simple_newick_summary(self, mocker, tmp_path, capsys):
        tree_path = tmp_path / "tree.tre"
        tree_path.write_text("(a:1,(b:2,c:3):4);")
        t = Treeness(Namespace(tree=str(tree_path), json=False))
        read_tree = mocker.patch.object(
            t,
            "read_tree_file_unmodified",
            side_effect=AssertionError("simple Newick should use summary path"),
        )

        t.run()

        read_tree.assert_not_called()
        out, _ = capsys.readouterr()
        assert out.strip() == "0.4"
