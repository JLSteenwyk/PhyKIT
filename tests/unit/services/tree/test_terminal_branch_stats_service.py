from argparse import Namespace
from io import StringIO
import subprocess
import sys

from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin
import pytest

from phykit.services.tree.terminal_branch_stats import TerminalBranchStats


def test_module_import_defers_json_helper_and_heavy_tree_dependencies():
    code = """
import importlib
import sys

module = importlib.import_module("phykit.services.tree.terminal_branch_stats")
assert callable(module.print_json)
assert hasattr(module.np, "__getattr__")
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.stats_summary" not in sys.modules
assert "typing" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "numpy" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


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

    def test_get_terminal_branch_lengths_can_skip_names(self, args):
        service = TerminalBranchStats(args)
        tree = _Tree([
            _Terminal("a", 1.0),
            _Terminal("b", None),
            _Terminal("c", 2.0),
        ])
        lengths, lengths_and_names = service.get_terminal_branch_lengths(
            tree, include_names=False
        )
        assert lengths == [1.0, 2.0]
        assert lengths_and_names == []

    def test_get_terminal_branch_lengths_uses_direct_traversal(self, monkeypatch, args):
        service = TerminalBranchStats(args)
        tree = Phylo.read(StringIO("((A:1,B:2):3,(C,D:4):5);"), "newick")

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("generic terminal traversal should not be used")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        lengths, lengths_and_names = service.get_terminal_branch_lengths(tree)

        assert lengths == [1.0, 2.0, 4.0]
        assert lengths_and_names == [[1.0, "A"], [2.0, "B"], [4.0, "D"]]

    def test_get_terminal_branch_lengths_without_names_uses_direct_traversal(
        self, monkeypatch, args
    ):
        service = TerminalBranchStats(args)
        tree = Phylo.read(StringIO("((A:1,B:2):3,(C,D:4):5);"), "newick")

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("generic terminal traversal should not be used")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        lengths, lengths_and_names = service.get_terminal_branch_lengths(
            tree,
            include_names=False,
        )

        assert lengths == [1.0, 2.0, 4.0]
        assert lengths_and_names == []

    def test_terminal_branch_lengths_direct_preserves_order_without_child_reversed(
        self, args
    ):
        class NoReversedList(list):
            def __reversed__(self):
                raise AssertionError("direct path should push children directly")

        class Clade:
            def __init__(self, branch_length=None, name=None, clades=None):
                self.branch_length = branch_length
                self.name = name
                self.clades = NoReversedList(clades or [])

        tree = type(
            "Tree",
            (),
            {
                "root": Clade(
                    clades=[
                        Clade(clades=[Clade(1.0, "A"), Clade(2.0, "B")]),
                        Clade(clades=[Clade(None, "C"), Clade(4.0, "D")]),
                    ],
                )
            },
        )()

        lengths, lengths_and_names = TerminalBranchStats._get_terminal_branch_lengths_direct(
            tree
        )

        assert lengths == [1.0, 2.0, 4.0]
        assert lengths_and_names == [[1.0, "A"], [2.0, "B"], [4.0, "D"]]

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

    def test_calculate_terminal_branch_stats_non_verbose_uses_array_path(
        self, mocker, args
    ):
        service = TerminalBranchStats(args)
        tree = Phylo.read(StringIO("((A:1,B:2):3,(C,D:4):5);"), "newick")
        mocker.patch.object(
            TerminalBranchStats,
            "get_terminal_branch_lengths",
            side_effect=AssertionError("non-verbose direct path should be used"),
        )
        mocked_stats = mocker.patch(
            "phykit.services.tree.terminal_branch_stats.calculate_summary_statistics_from_arr",
            return_value={"mean": 7 / 3},
        )

        lengths, stats, lengths_and_names = service.calculate_terminal_branch_stats(
            tree,
            include_names=False,
        )

        assert lengths == [1.0, 2.0, 4.0]
        assert stats == {"mean": 7 / 3}
        assert lengths_and_names == []
        stats_arg = mocked_stats.call_args.args[0]
        assert stats_arg.tolist() == [1.0, 2.0, 4.0]

    def test_terminal_branch_lengths_array_preserves_order_without_reversed(
        self, args
    ):
        class NoReversedList(list):
            def __reversed__(self):
                raise AssertionError("array path should push children directly")

        class Clade:
            def __init__(self, branch_length=None, clades=None):
                self.branch_length = branch_length
                self.clades = NoReversedList(clades or [])

        tree = type(
            "Tree",
            (),
            {
                "root": Clade(
                    clades=[
                        Clade(clades=[Clade(1.0), Clade(2.0)]),
                        Clade(clades=[Clade(None), Clade(4.0)]),
                    ],
                )
            },
        )()

        lengths = TerminalBranchStats._get_terminal_branch_lengths_array_direct(tree)

        assert lengths.tolist() == [1.0, 2.0, 4.0]

    def test_run_non_verbose_prints_summary(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", verbose=False, json=False)
        tree = _Tree([])
        service = TerminalBranchStats(args)
        mocker.patch.object(
            TerminalBranchStats,
            "read_tree_file_unmodified",
            return_value=tree,
        )
        mocked_stats = mocker.patch.object(
            TerminalBranchStats,
            "calculate_terminal_branch_stats",
            return_value=([1.0, 2.0], {"mean": 1.5}, [[1.0, "a"], [2.0, "b"]]),
        )
        mocked_print_summary = mocker.patch("phykit.services.tree.terminal_branch_stats.print_summary_statistics")

        service.run()
        mocked_print_summary.assert_called_once_with({"mean": 1.5})
        mocked_stats.assert_called_once_with(tree, include_names=False)

    def test_run_uses_unmodified_tree_read(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", verbose=False, json=False)
        tree = _Tree([])
        service = TerminalBranchStats(args)
        read_tree = mocker.patch.object(
            service,
            "read_tree_file_unmodified",
            return_value=tree,
        )
        mocker.patch.object(
            service,
            "calculate_terminal_branch_stats",
            return_value=([1.0], {"mean": 1.0}, []),
        )
        mocker.patch("phykit.services.tree.terminal_branch_stats.print_summary_statistics")

        service.run()

        read_tree.assert_called_once_with()

    def test_run_verbose_prints_lengths_and_names(self, mocker, capsys):
        args = Namespace(tree="/some/path/to/file.tre", verbose=True, json=False)
        tree = _Tree([])
        service = TerminalBranchStats(args)
        mocker.patch.object(
            TerminalBranchStats,
            "read_tree_file_unmodified",
            return_value=tree,
        )
        mocked_stats = mocker.patch.object(
            TerminalBranchStats,
            "calculate_terminal_branch_stats",
            return_value=([1.0, 2.0], {"mean": 1.5}, [[1.12345, "a"], [2.0, "b"]]),
        )

        service.run()
        out, _ = capsys.readouterr()
        assert out == "1.1235 a\n2.0 b\n"
        mocked_stats.assert_called_once_with(tree, include_names=True)

    def test_run_verbose_handles_broken_pipe(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", verbose=True, json=False)
        service = TerminalBranchStats(args)
        mocker.patch.object(
            TerminalBranchStats,
            "read_tree_file_unmodified",
            return_value=_Tree([]),
        )
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
        mocker.patch.object(
            TerminalBranchStats,
            "read_tree_file_unmodified",
            return_value=_Tree([]),
        )
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
        mocker.patch.object(
            TerminalBranchStats,
            "read_tree_file_unmodified",
            return_value=_Tree([]),
        )
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
