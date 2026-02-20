import pytest
from argparse import Namespace

from phykit.services.tree.treeness_over_rcv import TreenessOverRCV
import phykit.services.tree.treeness_over_rcv as tor_module


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
        # Mock the cached tree read method instead of Phylo.read
        mock_cached_read = mocker.patch("phykit.services.tree.base.Tree._cached_tree_read")
        mock_get_hash = mocker.patch("phykit.services.tree.base.Tree._get_file_hash", return_value="test_hash")

        t = TreenessOverRCV(args)
        t.read_tree_file()

        # Verify the cached read was called with the correct parameters
        mock_get_hash.assert_called_with(args.tree)
        mock_cached_read.assert_called_with(args.tree, "newick", "test_hash")

    def test_process_args_defaults_json_false(self):
        parsed = TreenessOverRCV(Namespace(tree="x.tre", alignment="x.fa")).process_args(
            Namespace(tree="x.tre", alignment="x.fa")
        )
        assert parsed["json_output"] is False

    def test_run_prints_tab_delimited(self, mocker, capsys):
        t = TreenessOverRCV(Namespace(tree="x.tre", alignment="x.fa", json=False))
        mocker.patch.object(TreenessOverRCV, "calculate_treeness", return_value=2.5)

        class FakeAlignment:
            def __init__(self, alignment_file_path):
                self.alignment_file_path = alignment_file_path

            def calculate_rcv(self):
                return 0.5

        mocker.patch.object(tor_module, "Alignment", FakeAlignment)
        t.run()
        out, _ = capsys.readouterr()
        assert out.strip() == "5.0\t2.5\t0.5"

    def test_run_json(self, mocker):
        t = TreenessOverRCV(Namespace(tree="x.tre", alignment="x.fa", json=True))
        mocker.patch.object(TreenessOverRCV, "calculate_treeness", return_value=1.23456)

        class FakeAlignment:
            def __init__(self, alignment_file_path):
                self.alignment_file_path = alignment_file_path

            def calculate_rcv(self):
                return 0.32109

        mocker.patch.object(tor_module, "Alignment", FakeAlignment)
        mocked_json = mocker.patch.object(tor_module, "print_json")
        t.run()
        mocked_json.assert_called_once_with(
            {"treeness_over_rcv": 3.8449, "treeness": 1.2346, "rcv": 0.3211}
        )
