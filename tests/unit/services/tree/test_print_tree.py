import pytest
from argparse import Namespace

from phykit.services.tree.print_tree import PrintTree
import phykit.services.tree.print_tree as print_tree_module


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

    def test_process_args_defaults_json_false(self):
        parsed = PrintTree(Namespace(tree="t.tre", remove=False)).process_args(
            Namespace(tree="t.tre", remove=False)
        )
        assert parsed["json_output"] is False

    def test_run_json_remove_branch_lengths(self, tree_simple, monkeypatch):
        captured = {}
        before_lengths = [node.branch_length for node in tree_simple.get_terminals()]
        t = PrintTree(Namespace(tree="x.tre", remove=True, json=True))
        monkeypatch.setattr(t, "read_tree_file", lambda: tree_simple)
        monkeypatch.setattr(print_tree_module, "print_json", lambda payload: captured.setdefault("payload", payload))

        t.run()

        assert captured["payload"]["remove_branch_lengths"] is True
        assert "0.00000" in captured["payload"]["tree_newick"]
        assert [node.branch_length for node in tree_simple.get_terminals()] == before_lengths

    def test_run_draw_ascii_and_broken_pipe(self, tree_simple, monkeypatch):
        t = PrintTree(Namespace(tree="x.tre", remove=False, json=False))
        monkeypatch.setattr(t, "read_tree_file", lambda: tree_simple)
        monkeypatch.setattr(print_tree_module.Phylo, "draw_ascii", lambda _tree: (_ for _ in ()).throw(BrokenPipeError()))

        # no exception should be raised
        t.run()
