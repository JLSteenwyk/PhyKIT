from argparse import Namespace
from io import StringIO
import subprocess
import sys

import pytest
from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin

import phykit.services.tree.spurious_sequence as spurious_sequence_module
from phykit.services.tree.spurious_sequence import SpuriousSequence


def test_module_import_does_not_import_statistics_or_json():
    code = """
import sys
import phykit.services.tree.spurious_sequence as module
assert callable(module.print_json)
assert hasattr(module.np, "__getattr__")
assert "typing" not in sys.modules
assert "statistics" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "numpy" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


@pytest.fixture
def args():
    return Namespace(tree="/some/path/to/file.tre", factor=None)


class _Terminal:
    def __init__(self, name, branch_length):
        self.name = name
        self.branch_length = branch_length


class _Tree:
    def __init__(self, terminals):
        self._terminals = terminals

    def get_terminals(self):
        return self._terminals


class TestSpuriousSequence:
    def test_init_sets_expected_attrs(self, args):
        service = SpuriousSequence(args)
        assert service.tree_file_path == args.tree
        assert service.factor == 20
        assert service.json_output is False

    def test_get_branch_lengths_and_names_uses_terminal_branches_only(self, args):
        service = SpuriousSequence(args)
        tree = _Tree(
            [
                _Terminal("a", 1.0),
                _Terminal("b", None),
                _Terminal("c", 3.0),
            ]
        )
        branch_lengths, name_map = service.get_branch_lengths_and_their_names(tree)
        assert branch_lengths == [1.0, 3.0]
        assert name_map == {"a": 1.0, "c": 3.0}

    def test_get_branch_lengths_fast_path_does_not_call_get_terminals(
        self, args, monkeypatch
    ):
        service = SpuriousSequence(args)
        tree = Phylo.read(StringIO("((a:1,b:2):3,c:4);"), "newick")

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("standard tree should use direct traversal")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        branch_lengths, name_map = service.get_branch_lengths_and_their_names(tree)

        assert branch_lengths == [1.0, 2.0, 4.0]
        assert name_map == {"a": 1.0, "b": 2.0, "c": 4.0}

    def test_get_branch_lengths_standard_tree_collects_in_one_pass(
        self, args, monkeypatch
    ):
        service = SpuriousSequence(args)
        tree = Phylo.read(StringIO("((a:1,b:2):3,c:4);"), "newick")

        def fail_terminal_list(_tree):
            raise AssertionError(
                "standard tree should collect branch lengths directly"
            )

        monkeypatch.setattr(
            SpuriousSequence,
            "_iter_terminal_clades",
            staticmethod(fail_terminal_list),
        )

        branch_lengths, name_map = service.get_branch_lengths_and_their_names(tree)

        assert branch_lengths == [1.0, 2.0, 4.0]
        assert name_map == {"a": 1.0, "b": 2.0, "c": 4.0}

    def test_iter_terminal_clades_preserves_order_with_mixed_child_counts(
        self, monkeypatch
    ):
        tree = Phylo.read(
            StringIO("(a:1,(b:2,c:3):4,(d:5,e:6,f:7):8);"),
            "newick",
        )

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("standard tree should use direct traversal")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        terminals = SpuriousSequence._iter_terminal_clades(tree)

        assert [terminal.name for terminal in terminals] == [
            "a",
            "b",
            "c",
            "d",
            "e",
            "f",
        ]

    def test_identify_spurious_sequence(self, args):
        service = SpuriousSequence(args)
        tree = _Tree([_Terminal("a", 1.0), _Terminal("b", 3.0), _Terminal("c", 5.0)])
        name_map, threshold, median = service.identify_spurious_sequence(tree, factor=2)
        assert median == 3.0
        assert threshold == 6.0
        assert name_map["c"] == 5.0

    def test_median_branch_length_matches_statistics_median(self):
        assert SpuriousSequence._median_branch_length([3.0, 1.0, 5.0]) == 3.0
        assert SpuriousSequence._median_branch_length([4.0, 1.0, 3.0, 2.0]) == 2.5

    def test_median_branch_length_uses_partition_for_large_inputs(self, monkeypatch):
        calls = []

        class _FakeArray:
            def __init__(self, values):
                self.values = sorted(values)

            def partition(self, kth):
                calls.append(kth)

            def __getitem__(self, index):
                return self.values[index]

        class _FakeNumpy:
            def asarray(self, values, dtype=float):
                assert dtype is float
                return _FakeArray(values)

        monkeypatch.setattr(spurious_sequence_module, "_MEDIAN_NUMPY_THRESHOLD", 4)
        monkeypatch.setattr(spurious_sequence_module, "np", _FakeNumpy())

        assert SpuriousSequence._median_branch_length([4.0, 1.0, 3.0, 2.0]) == 2.5
        assert calls == [(1, 2)]

    def test_median_branch_length_default_threshold_uses_partition(self, monkeypatch):
        calls = []

        class _FakeArray:
            def __init__(self, values):
                self.values = sorted(values)

            def partition(self, kth):
                calls.append(kth)

            def __getitem__(self, index):
                return self.values[index]

        class _FakeNumpy:
            def asarray(self, values, dtype=float):
                assert dtype is float
                return _FakeArray(values)

        monkeypatch.setattr(spurious_sequence_module, "np", _FakeNumpy())

        values = list(range(spurious_sequence_module._MEDIAN_NUMPY_THRESHOLD))
        assert SpuriousSequence._median_branch_length(values) == 511.5
        assert calls == [(511, 512)]

    def test_run_prints_none_when_no_spurious_sequences(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", factor=20, json=False)
        service = SpuriousSequence(args)
        mocker.patch.object(
            SpuriousSequence,
            "read_tree_file_unmodified",
            return_value=_Tree([]),
        )
        mocker.patch.object(
            SpuriousSequence,
            "identify_spurious_sequence",
            return_value=({"a": 1.0, "b": 2.0}, 100.0, 1.5),
        )
        mocked_print = mocker.patch("builtins.print")
        service.run()
        mocked_print.assert_called_once_with("None")

    def test_run_prints_spurious_sequences(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", factor=2, json=False)
        service = SpuriousSequence(args)
        mocker.patch.object(
            SpuriousSequence,
            "read_tree_file_unmodified",
            return_value=_Tree([]),
        )
        mocker.patch.object(
            SpuriousSequence,
            "identify_spurious_sequence",
            return_value=({"a": 10.0, "b": 8.0, "c": 1.0}, 5.0, 2.5),
        )
        mocked_print = mocker.patch("builtins.print")
        service.run()
        mocked_print.assert_called_once_with("a\t10.0\t5.0\t2.5\nb\t8.0\t5.0\t2.5")

    def test_run_json_payload(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", factor=2, json=True)
        service = SpuriousSequence(args)
        mocker.patch.object(
            SpuriousSequence,
            "read_tree_file_unmodified",
            return_value=_Tree([]),
        )
        mocker.patch.object(
            SpuriousSequence,
            "identify_spurious_sequence",
            return_value=({"a": 10.12345, "b": 1.0}, 5.0, 2.5),
        )
        mocked_json = mocker.patch("phykit.services.tree.spurious_sequence.print_json")
        service.run()
        payload = mocked_json.call_args.args[0]
        assert payload["rows"] == payload["spurious_sequences"]
        assert payload["rows"][0] == {
            "taxon": "a",
            "branch_length": 10.1235,
            "threshold": 5.0,
            "median": 2.5,
        }

    def test_run_uses_unmodified_tree_read(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", factor=2, json=False)
        tree = _Tree([])
        service = SpuriousSequence(args)
        read_tree = mocker.patch.object(
            service,
            "read_tree_file_unmodified",
            return_value=tree,
        )
        mocker.patch.object(
            service,
            "identify_spurious_sequence",
            return_value=({"a": 1.0}, 10.0, 1.0),
        )
        mocked_print = mocker.patch("builtins.print")

        service.run()

        read_tree.assert_called_once_with()
        mocked_print.assert_called_once_with("None")
