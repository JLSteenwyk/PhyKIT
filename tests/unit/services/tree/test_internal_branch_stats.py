import pytest
from argparse import Namespace
from math import isclose

from phykit.services.tree.internal_branch_stats import InternalBranchStats
import phykit.services.tree.internal_branch_stats as internal_branch_stats_module


@pytest.fixture
def args():
    kwargs = dict(tree="/some/path/to/file.tre", verbose=None)
    return Namespace(**kwargs)


class TestInternalBranchStats(object):
    def test_init_sets_tree_file_path(self, args):
        t = InternalBranchStats(args)
        assert t.tree_file_path == args.tree
        assert t.output_file_path is None

    def test_read_file_reads_tree_file_path(self, mocker, args):
        # Mock the cached tree read method instead of Phylo.read
        mock_cached_read = mocker.patch("phykit.services.tree.base.Tree._cached_tree_read")
        mock_get_hash = mocker.patch("phykit.services.tree.base.Tree._get_file_hash", return_value="test_hash")

        t = InternalBranchStats(args)
        t.read_tree_file()

        # Verify the cached read was called with the correct parameters
        mock_get_hash.assert_called_with(args.tree)
        mock_cached_read.assert_called_with(args.tree, "newick", "test_hash")

    def test_calculate_internal_branch_stats(self, tree_simple, args):
        t = InternalBranchStats(args)
        stats, _ = t.calculate_internal_branch_stats(tree_simple)
        assert isinstance(stats["mean"], float)
        assert isinstance(stats["median"], float)
        assert isinstance(stats["twenty_fifth"], float)
        assert isinstance(stats["seventy_fifth"], float)
        assert isinstance(stats["standard_deviation"], float)
        assert isinstance(stats["variance"], float)
        assert isinstance(stats["minimum"], float)
        assert isinstance(stats["maximum"], float)
        assert isclose(stats["mean"], 6.987232, rel_tol=0.001)
        assert isclose(stats["median"], 3.87382, rel_tol=0.001)
        assert isclose(stats["twenty_fifth"], 2.0946, rel_tol=0.001)
        assert isclose(stats["seventy_fifth"], 7.52973, rel_tol=0.001)
        assert isclose(stats["standard_deviation"], 8.011401268758792, rel_tol=0.001)
        assert isclose(stats["variance"], 64.18255028907, rel_tol=0.001)
        assert isclose(stats["minimum"], 0.846, rel_tol=0.001)
        assert isclose(stats["maximum"], 20.59201, rel_tol=0.001)

    def test_process_args_defaults_json_false(self):
        parsed = InternalBranchStats(Namespace(tree="x.tre", verbose=False)).process_args(
            Namespace(tree="x.tre", verbose=False)
        )
        assert parsed["json_output"] is False

    def test_get_internal_branch_lengths(self, tree_simple, args):
        t = InternalBranchStats(args)
        lengths, lengths_and_names = t.get_internal_branch_lengths(tree_simple)
        assert len(lengths) > 0
        assert len(lengths) == len(lengths_and_names)
        assert isinstance(lengths_and_names[0][1], list)

    def test_calculate_internal_branch_stats_without_lengths_exits(self, capsys):
        from Bio import Phylo
        from io import StringIO

        t = InternalBranchStats(Namespace(tree="x.tre", verbose=False))
        tree = Phylo.read(StringIO("(A,B,C);"), "newick")
        with pytest.raises(SystemExit) as exc:
            t.calculate_internal_branch_stats(tree)
        assert exc.value.code == 2
        out, _ = capsys.readouterr()
        assert "requires a phylogeny with branch lengths" in out

    def test_run_verbose_json(self, mocker):
        t = InternalBranchStats(Namespace(tree="x.tre", verbose=True, json=True))
        mocker.patch.object(InternalBranchStats, "read_tree_file", return_value=object())
        mocker.patch.object(
            InternalBranchStats,
            "calculate_internal_branch_stats",
            return_value=({"mean": 1.0}, [(1.23456, ["A", "B"])]),
        )
        mocked_json = mocker.patch.object(internal_branch_stats_module, "print_json")
        t.run()
        payload = mocked_json.call_args.args[0]
        assert payload["verbose"] is True
        assert payload["rows"] == payload["branches"]
        assert payload["rows"][0] == {"length": 1.2346, "terminals": ["A", "B"]}

    def test_run_nonverbose_json(self, mocker):
        t = InternalBranchStats(Namespace(tree="x.tre", verbose=False, json=True))
        mocker.patch.object(InternalBranchStats, "read_tree_file", return_value=object())
        mocker.patch.object(
            InternalBranchStats,
            "calculate_internal_branch_stats",
            return_value=({"mean": 1.0}, [(1.2, ["A"])]),
        )
        mocked_json = mocker.patch.object(internal_branch_stats_module, "print_json")
        t.run()
        mocked_json.assert_called_once_with({"verbose": False, "summary": {"mean": 1.0}})

    def test_run_verbose_print_handles_broken_pipe(self, mocker):
        t = InternalBranchStats(Namespace(tree="x.tre", verbose=True, json=False))
        mocker.patch.object(InternalBranchStats, "read_tree_file", return_value=object())
        mocker.patch.object(
            InternalBranchStats,
            "calculate_internal_branch_stats",
            return_value=({"mean": 1.0}, [(1.2, ["A", "B"])]),
        )
        mocker.patch("builtins.print", side_effect=BrokenPipeError)
        t.run()

    def test_run_nonverbose_prints_summary(self, mocker):
        t = InternalBranchStats(Namespace(tree="x.tre", verbose=False, json=False))
        mocker.patch.object(InternalBranchStats, "read_tree_file", return_value=object())
        mocker.patch.object(
            InternalBranchStats,
            "calculate_internal_branch_stats",
            return_value=({"mean": 1.0}, [(1.2, ["A", "B"])]),
        )
        mocked_summary = mocker.patch.object(internal_branch_stats_module, "print_summary_statistics")
        t.run()
        mocked_summary.assert_called_once_with({"mean": 1.0})
