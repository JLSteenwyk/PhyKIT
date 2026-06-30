import pytest
import subprocess
import sys
from argparse import Namespace
from math import isclose
from io import StringIO

from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin

from phykit.services.tree.internal_branch_stats import InternalBranchStats
import phykit.services.tree.internal_branch_stats as internal_branch_stats_module


def test_module_import_defers_json_helper_and_heavy_tree_dependencies():
    code = """
import importlib
import sys

module = importlib.import_module("phykit.services.tree.internal_branch_stats")
assert callable(module.print_json)
assert hasattr(module.np, "__getattr__")
assert "json" not in sys.modules
assert "typing" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.stats_summary" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "numpy" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


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

    def test_get_internal_branch_lengths_matches_terminal_scan(self, tree_simple, args):
        t = InternalBranchStats(args)
        lengths, lengths_and_names = t.get_internal_branch_lengths(tree_simple)

        expected_lengths = []
        expected_rows = []
        for internal_branch in tree_simple.get_nonterminals():
            if internal_branch.branch_length is not None:
                expected_lengths.append(internal_branch.branch_length)
                expected_rows.append(
                    (
                        internal_branch.branch_length,
                        [term.name for term in internal_branch.get_terminals()],
                    )
                )

        assert lengths == expected_lengths
        assert lengths_and_names == expected_rows

    def test_get_internal_branch_lengths_can_skip_terminal_names(self, tree_simple, args):
        t = InternalBranchStats(args)
        lengths_with_names, _ = t.get_internal_branch_lengths(tree_simple)
        lengths, lengths_and_names = t.get_internal_branch_lengths(
            tree_simple, include_names=False
        )

        assert lengths == lengths_with_names
        assert lengths_and_names == []

    def test_get_internal_branch_lengths_uses_direct_traversal(self, monkeypatch, args):
        t = InternalBranchStats(args)
        tree = Phylo.read(StringIO("((A:1,B:2):3,(C:4,D:5):6):7;"), "newick")

        def fail_generic_traversal(*_args, **_kwargs):
            raise AssertionError("generic tree traversal should not be used")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_generic_traversal)
        monkeypatch.setattr(TreeMixin, "get_nonterminals", fail_generic_traversal)

        lengths, lengths_and_names = t.get_internal_branch_lengths(tree)

        assert lengths == [7.0, 3.0, 6.0]
        assert lengths_and_names == [
            (7.0, ["A", "B", "C", "D"]),
            (3.0, ["A", "B"]),
            (6.0, ["C", "D"]),
        ]

    def test_get_internal_branch_lengths_matches_generic_reference_with_direct_path(
        self, monkeypatch, args
    ):
        t = InternalBranchStats(args)
        tree = Phylo.read(
            StringIO("((A:1,B:2):3,(C:4,(D:5,E:6):7):8):9;"),
            "newick",
        )
        expected_lengths = []
        expected_rows = []
        for internal_branch in tree.get_nonterminals():
            if internal_branch.branch_length is not None:
                expected_lengths.append(internal_branch.branch_length)
                expected_rows.append(
                    (
                        internal_branch.branch_length,
                        [term.name for term in internal_branch.get_terminals()],
                    )
                )

        def fail_generic_traversal(*_args, **_kwargs):
            raise AssertionError("standard tree helper should use direct traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_generic_traversal)
        monkeypatch.setattr(TreeMixin, "get_nonterminals", fail_generic_traversal)

        lengths, rows = t.get_internal_branch_lengths(tree)

        assert lengths == expected_lengths
        assert rows == expected_rows

    def test_get_internal_branch_lengths_without_names_uses_direct_traversal(
        self, monkeypatch, args
    ):
        t = InternalBranchStats(args)
        tree = Phylo.read(StringIO("((A:1,B:2):3,(C:4,D:5):6):7;"), "newick")

        def fail_generic_traversal(*_args, **_kwargs):
            raise AssertionError("generic tree traversal should not be used")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_generic_traversal)
        monkeypatch.setattr(TreeMixin, "get_nonterminals", fail_generic_traversal)

        lengths, lengths_and_names = t.get_internal_branch_lengths(
            tree, include_names=False
        )

        assert lengths == [7.0, 3.0, 6.0]
        assert lengths_and_names == []

    def test_get_internal_branch_lengths_without_names_handles_mixed_child_counts(
        self, monkeypatch, args
    ):
        t = InternalBranchStats(args)
        tree = Phylo.read(
            StringIO("(A:1,(B:2,C:3,D:4):5,(E:6,F:7):8):9;"),
            "newick",
        )

        def fail_generic_traversal(*_args, **_kwargs):
            raise AssertionError("standard tree helper should use direct traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_generic_traversal)
        monkeypatch.setattr(TreeMixin, "get_nonterminals", fail_generic_traversal)

        lengths, lengths_and_names = t.get_internal_branch_lengths(
            tree, include_names=False
        )

        assert lengths == [9.0, 5.0, 8.0]
        assert lengths_and_names == []

    def test_internal_branch_lengths_verbose_preserves_order_without_reversed(self):
        class NoReversedList(list):
            def __reversed__(self):
                raise AssertionError("verbose path should push children directly")

        class Clade:
            def __init__(self, name=None, branch_length=None, clades=None):
                self.name = name
                self.branch_length = branch_length
                self.clades = NoReversedList(clades or [])

        tree = type(
            "Tree",
            (),
            {
                "root": Clade(
                    branch_length=7.0,
                    clades=[
                        Clade(
                            branch_length=3.0,
                            clades=[Clade("A", 1.0), Clade("B", 2.0)],
                        ),
                        Clade(
                            branch_length=6.0,
                            clades=[Clade("C", 4.0), Clade("D", 5.0)],
                        ),
                    ],
                )
            },
        )()

        lengths, rows = InternalBranchStats._get_internal_branch_lengths_direct(
            tree,
            include_names=True,
        )

        assert lengths == [7.0, 3.0, 6.0]
        assert rows == [
            (7.0, ["A", "B", "C", "D"]),
            (3.0, ["A", "B"]),
            (6.0, ["C", "D"]),
        ]

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

    def test_calculate_internal_branch_stats_nonverbose_uses_array_path(
        self, mocker, args
    ):
        t = InternalBranchStats(args)
        tree = Phylo.read(StringIO("((A:1,B:2):3,(C:4,D:5):6):7;"), "newick")
        mocker.patch.object(
            InternalBranchStats,
            "get_internal_branch_lengths",
            side_effect=AssertionError("non-verbose direct path should be used"),
        )
        mocked_stats = mocker.patch(
            "phykit.services.tree.internal_branch_stats.calculate_summary_statistics_from_arr",
            return_value={"mean": 16 / 3},
        )

        stats, lengths_and_names = t.calculate_internal_branch_stats(
            tree,
            include_names=False,
        )

        assert stats == {"mean": 16 / 3}
        assert lengths_and_names == []
        stats_arg = mocked_stats.call_args.args[0]
        assert stats_arg.tolist() == [7.0, 3.0, 6.0]

    def test_internal_branch_lengths_array_preserves_order_without_reversed(self):
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
                    7.0,
                    [
                        Clade(3.0, [Clade(1.0), Clade(2.0)]),
                        Clade(6.0, [Clade(4.0), Clade(5.0)]),
                    ],
                )
            },
        )()

        lengths = InternalBranchStats._get_internal_branch_lengths_array_direct(tree)

        assert lengths.tolist() == [7.0, 3.0, 6.0]

    def test_run_verbose_json(self, mocker):
        tree = object()
        t = InternalBranchStats(Namespace(tree="x.tre", verbose=True, json=True))
        mocker.patch.object(
            InternalBranchStats,
            "read_tree_file_unmodified",
            return_value=tree,
        )
        mocked_stats = mocker.patch.object(
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
        mocked_stats.assert_called_once_with(tree, include_names=True)

    def test_run_nonverbose_json(self, mocker):
        tree = object()
        t = InternalBranchStats(Namespace(tree="x.tre", verbose=False, json=True))
        mocker.patch.object(
            InternalBranchStats,
            "read_tree_file_unmodified",
            return_value=tree,
        )
        mocked_stats = mocker.patch.object(
            InternalBranchStats,
            "calculate_internal_branch_stats",
            return_value=({"mean": 1.0}, [(1.2, ["A"])]),
        )
        mocked_json = mocker.patch.object(internal_branch_stats_module, "print_json")
        t.run()
        mocked_json.assert_called_once_with({"verbose": False, "summary": {"mean": 1.0}})
        mocked_stats.assert_called_once_with(tree, include_names=False)

    def test_run_uses_unmodified_tree_read(self, mocker):
        tree = object()
        t = InternalBranchStats(Namespace(tree="x.tre", verbose=False, json=False))
        read_tree = mocker.patch.object(
            t,
            "read_tree_file_unmodified",
            return_value=tree,
        )
        mocker.patch.object(
            t,
            "calculate_internal_branch_stats",
            return_value=({"mean": 1.0}, []),
        )
        mocker.patch.object(internal_branch_stats_module, "print_summary_statistics")

        t.run()

        read_tree.assert_called_once_with()

    def test_run_verbose_prints_lengths_and_names(self, mocker, capsys):
        tree = object()
        t = InternalBranchStats(Namespace(tree="x.tre", verbose=True, json=False))
        mocker.patch.object(
            InternalBranchStats,
            "read_tree_file_unmodified",
            return_value=tree,
        )
        mocked_stats = mocker.patch.object(
            InternalBranchStats,
            "calculate_internal_branch_stats",
            return_value=({"mean": 1.0}, [(1.23456, ["A", "B"]), (2.0, ["C"])]),
        )

        t.run()

        out, _ = capsys.readouterr()
        assert out == "1.2346 A;B\n2.0 C\n"
        mocked_stats.assert_called_once_with(tree, include_names=True)

    def test_run_verbose_print_handles_broken_pipe(self, mocker):
        t = InternalBranchStats(Namespace(tree="x.tre", verbose=True, json=False))
        mocker.patch.object(
            InternalBranchStats,
            "read_tree_file_unmodified",
            return_value=object(),
        )
        mocker.patch.object(
            InternalBranchStats,
            "calculate_internal_branch_stats",
            return_value=({"mean": 1.0}, [(1.2, ["A", "B"])]),
        )
        mocker.patch("builtins.print", side_effect=BrokenPipeError)
        t.run()

    def test_run_nonverbose_prints_summary(self, mocker):
        t = InternalBranchStats(Namespace(tree="x.tre", verbose=False, json=False))
        mocker.patch.object(
            InternalBranchStats,
            "read_tree_file_unmodified",
            return_value=object(),
        )
        mocker.patch.object(
            InternalBranchStats,
            "calculate_internal_branch_stats",
            return_value=({"mean": 1.0}, [(1.2, ["A", "B"])]),
        )
        mocked_summary = mocker.patch.object(internal_branch_stats_module, "print_summary_statistics")
        t.run()
        mocked_summary.assert_called_once_with({"mean": 1.0})
