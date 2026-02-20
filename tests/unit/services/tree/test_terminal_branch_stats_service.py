from argparse import Namespace

import pytest

from phykit.services.tree.terminal_branch_stats import TerminalBranchStats


@pytest.fixture
def args():
    return Namespace(tree="/some/path/to/file.tre", verbose=False)


class _Terminal:
    def __init__(self, name, branch_length):
        self.name = name
        self.branch_length = branch_length


class _Tree:
    def __init__(self, terminals):
        self._terminals = terminals

    def get_terminals(self):
        return self._terminals


class TestTerminalBranchStats:
    def test_init_sets_expected_attrs(self, args):
        service = TerminalBranchStats(args)
        assert service.tree_file_path == args.tree
        assert service.verbose is False
        assert service.json_output is False

    def test_get_terminal_branch_lengths_skips_none(self, args):
        service = TerminalBranchStats(args)
        tree = _Tree([
            _Terminal("a", 1.0),
            _Terminal("b", None),
            _Terminal("c", 2.0),
        ])
        lengths, lengths_and_names = service.get_terminal_branch_lengths(tree)
        assert lengths == [1.0, 2.0]
        assert lengths_and_names == [[1.0, "a"], [2.0, "c"]]

    def test_check_tree_has_branch_lengths_exits_on_empty(self, args):
        service = TerminalBranchStats(args)
        with pytest.raises(SystemExit) as excinfo:
            service.check_tree_has_branch_lengths([])
        assert excinfo.value.code == 2

    def test_calculate_terminal_branch_stats(self, mocker, args):
        service = TerminalBranchStats(args)
        tree = _Tree([_Terminal("a", 1.0), _Terminal("b", 2.0)])
        mocker.patch(
            "phykit.services.tree.terminal_branch_stats.calculate_summary_statistics_from_arr",
            return_value={"mean": 1.5},
        )
        lengths, stats, lengths_and_names = service.calculate_terminal_branch_stats(tree)
        assert lengths == [1.0, 2.0]
        assert stats == {"mean": 1.5}
        assert lengths_and_names == [[1.0, "a"], [2.0, "b"]]

    def test_run_non_verbose_prints_summary(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", verbose=False, json=False)
        service = TerminalBranchStats(args)
        mocker.patch.object(TerminalBranchStats, "read_tree_file", return_value=_Tree([]))
        mocker.patch.object(
            TerminalBranchStats,
            "calculate_terminal_branch_stats",
            return_value=([1.0, 2.0], {"mean": 1.5}, [[1.0, "a"], [2.0, "b"]]),
        )
        mocked_print_summary = mocker.patch("phykit.services.tree.terminal_branch_stats.print_summary_statistics")

        service.run()
        mocked_print_summary.assert_called_once_with({"mean": 1.5})

    def test_run_verbose_prints_lengths_and_names(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", verbose=True, json=False)
        service = TerminalBranchStats(args)
        mocker.patch.object(TerminalBranchStats, "read_tree_file", return_value=_Tree([]))
        mocker.patch.object(
            TerminalBranchStats,
            "calculate_terminal_branch_stats",
            return_value=([1.0, 2.0], {"mean": 1.5}, [[1.12345, "a"], [2.0, "b"]]),
        )
        mocked_print = mocker.patch("builtins.print")

        service.run()
        mocked_print.assert_any_call(1.1235, "a")
        mocked_print.assert_any_call(2.0, "b")

    def test_run_verbose_handles_broken_pipe(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", verbose=True, json=False)
        service = TerminalBranchStats(args)
        mocker.patch.object(TerminalBranchStats, "read_tree_file", return_value=_Tree([]))
        mocker.patch.object(
            TerminalBranchStats,
            "calculate_terminal_branch_stats",
            return_value=([1.0], {"mean": 1.0}, [[1.0, "a"]]),
        )
        mocker.patch("builtins.print", side_effect=BrokenPipeError)

        service.run()

    def test_run_json_summary(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", verbose=False, json=True)
        service = TerminalBranchStats(args)
        mocker.patch.object(TerminalBranchStats, "read_tree_file", return_value=_Tree([]))
        mocker.patch.object(
            TerminalBranchStats,
            "calculate_terminal_branch_stats",
            return_value=([1.0, 2.0], {"mean": 1.5}, [[1.0, "a"], [2.0, "b"]]),
        )
        mocked_json = mocker.patch("phykit.services.tree.terminal_branch_stats.print_json")

        service.run()
        payload = mocked_json.call_args.args[0]
        assert payload == {"verbose": False, "summary": {"mean": 1.5}}

    def test_run_json_verbose(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", verbose=True, json=True)
        service = TerminalBranchStats(args)
        mocker.patch.object(TerminalBranchStats, "read_tree_file", return_value=_Tree([]))
        mocker.patch.object(
            TerminalBranchStats,
            "calculate_terminal_branch_stats",
            return_value=([1.0, 2.0], {"mean": 1.5}, [[1.12345, "a"], [2.0, "b"]]),
        )
        mocked_json = mocker.patch("phykit.services.tree.terminal_branch_stats.print_json")

        service.run()
        payload = mocked_json.call_args.args[0]
        assert payload["verbose"] is True
        assert payload["rows"] == [{"length": 1.1235, "taxon": "a"}, {"length": 2.0, "taxon": "b"}]
        assert payload["tips"] == payload["rows"]
