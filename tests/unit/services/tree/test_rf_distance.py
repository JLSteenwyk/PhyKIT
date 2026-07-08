import copy
import pickle as stdlib_pickle
import pytest
import subprocess
import sys
from argparse import Namespace
from concurrent.futures import Future
from io import StringIO

from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin

import phykit.services.tree.rf_distance as rf_distance_module
from phykit.services.tree.rf_distance import RobinsonFouldsDistance


def test_module_import_does_not_import_biophylo_or_numpy():
    code = """
import sys
import phykit.services.tree.rf_distance

assert "typing" not in sys.modules
assert "pickle" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "numpy" not in sys.modules
assert "concurrent.futures" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_lazy_pickle_caches_resolved_batch_helpers(monkeypatch):
    calls = []

    def cached_dumps(*args, **kwargs):
        calls.append(("dumps", args, kwargs))
        return args[0]

    def cached_loads(*args, **kwargs):
        calls.append(("loads", args, kwargs))
        return args[0]

    def uncached_dumps(*_args, **_kwargs):
        return "uncached-dumps"

    def uncached_loads(*_args, **_kwargs):
        return "uncached-loads"

    lazy_pickle = rf_distance_module._LazyPickle()
    monkeypatch.setattr(stdlib_pickle, "dumps", cached_dumps)
    monkeypatch.setattr(stdlib_pickle, "loads", cached_loads)

    assert lazy_pickle.loads(lazy_pickle.dumps("batch")) == "batch"

    monkeypatch.setattr(stdlib_pickle, "dumps", uncached_dumps)
    monkeypatch.setattr(stdlib_pickle, "loads", uncached_loads)

    assert lazy_pickle.loads(lazy_pickle.dumps("batch2")) == "batch2"
    assert lazy_pickle.__dict__["dumps"] is cached_dumps
    assert lazy_pickle.__dict__["loads"] is cached_loads
    assert calls == [
        ("dumps", ("batch",), {}),
        ("loads", ("batch",), {}),
        ("dumps", ("batch2",), {}),
        ("loads", ("batch2",), {}),
    ]


@pytest.fixture
def args():
    kwargs = dict(
        tree_zero="/some/path/to/file_zero.tre",
        tree_one="/some/path/to/file_one.tre",
    )
    return Namespace(**kwargs)


class TestRobinsonFouldsDistance(object):
    def test_init_sets_tree_file_path(self, args):
        rf = RobinsonFouldsDistance(args)
        assert rf.tree_file_path == args.tree_zero
        assert rf.tree1_file_path == args.tree_one
        assert rf.output_file_path is None

    def test_read_file_reads_tree_file_path(self, mocker, args):
        # Mock the cached tree read method instead of Phylo.read
        mock_cached_read = mocker.patch("phykit.services.tree.base.Tree._cached_tree_read")
        mock_get_hash = mocker.patch("phykit.services.tree.base.Tree._get_file_hash", return_value="test_hash")

        rf = RobinsonFouldsDistance(args)
        rf.read_tree_file()

        # Verify the cached read was called with the correct parameters
        mock_get_hash.assert_called_with(args.tree_zero)
        mock_cached_read.assert_called_with(args.tree_zero, "newick", "test_hash")

    def test_calculate_patristic_distances(self, tree_simple, tree_simple_other, args):
        rf = RobinsonFouldsDistance(args)
        plain_rf, normalized_rf = rf.calculate_robinson_foulds_distance(
            tree_simple, tree_simple_other
        )
        assert isinstance(plain_rf, int)
        assert isinstance(normalized_rf, float)
        assert plain_rf == 8
        assert normalized_rf == 0.8

    def test_calculate_robinson_foulds_distance_identical_trees(self, tree_simple, args):
        rf = RobinsonFouldsDistance(args)
        tree_clone = copy.deepcopy(tree_simple)

        plain_rf, normalized_rf = rf.calculate_robinson_foulds_distance(tree_simple, tree_clone)

        assert plain_rf == 0
        assert normalized_rf == 0

    def test_calculate_robinson_foulds_distance_same_object_skips_bipartitions(
        self, mocker, args
    ):
        rf = RobinsonFouldsDistance(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        mocker.patch.object(
            rf,
            "_get_all_bipartition_id_sets_direct",
            side_effect=AssertionError("same-object RF should skip bipartitions"),
        )

        assert rf.calculate_robinson_foulds_distance(tree, tree) == (0, 0.0)

    def test_calculate_robinson_foulds_distance_same_three_tip_object_preserves_denominator(
        self, args
    ):
        rf = RobinsonFouldsDistance(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,C:1);"), "newick")

        with pytest.raises(ZeroDivisionError):
            rf.calculate_robinson_foulds_distance(tree, tree)

    def test_get_all_bipartitions_matches_legacy_tree_comparison(
        self, tree_simple, tree_simple_other, args
    ):
        rf = RobinsonFouldsDistance(args)

        legacy_rf = 0
        legacy_rf = rf.compare_trees_optimized(legacy_rf, tree_simple, tree_simple_other)
        legacy_rf = rf.compare_trees_optimized(legacy_rf, tree_simple_other, tree_simple)

        split_rf = len(
            rf.get_all_bipartitions(tree_simple)
            ^ rf.get_all_bipartitions(tree_simple_other)
        )

        assert split_rf == legacy_rf

    def test_get_all_bipartitions_uses_direct_postorder(
        self, monkeypatch, args
    ):
        rf = RobinsonFouldsDistance(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("generic postorder traversal should not be used")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        bipartitions = rf.get_all_bipartitions(tree)

        assert bipartitions == {
            frozenset({"A", "B"}),
            frozenset({"C", "D"}),
        }

    def test_bipartition_id_sets_match_public_bipartitions(self, args):
        rf = RobinsonFouldsDistance(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")

        id_sets, tip_index = rf._get_all_bipartition_id_sets_direct(tree)
        index_to_tip = {idx: tip for tip, idx in tip_index.items()}
        named_sets = {
            frozenset(index_to_tip[idx] for idx in split)
            for split in id_sets
        }

        assert named_sets == rf.get_all_bipartitions(tree)

    def test_direct_bipartitions_handle_mixed_child_counts(self, args):
        rf = RobinsonFouldsDistance(args)
        tree = Phylo.read(
            StringIO("(((A:1):1,(B:1,C:1,D:1):1):1,E:1);"),
            "newick",
        )

        assert rf._get_all_bipartitions_direct(tree) == {
            frozenset({"A"}),
            frozenset({"B", "C", "D"}),
            frozenset({"A", "B", "C", "D"}),
        }

    def test_direct_bipartition_id_sets_handle_mixed_child_counts(self, args):
        rf = RobinsonFouldsDistance(args)
        tree = Phylo.read(
            StringIO("(((A:1):1,(B:1,C:1,D:1):1):1,E:1);"),
            "newick",
        )

        id_sets, tip_index = rf._get_all_bipartition_id_sets_direct(tree)
        index_to_tip = {idx: tip for tip, idx in tip_index.items()}
        named_sets = {
            frozenset(index_to_tip[idx] for idx in split)
            for split in id_sets
        }

        assert named_sets == {
            frozenset({"A"}),
            frozenset({"B", "C", "D"}),
            frozenset({"A", "B", "C", "D"}),
        }

    def test_direct_bipartition_helpers_fallback_on_malformed_clade(self):
        class MalformedChild:
            name = "A"

        class Root:
            clades = [MalformedChild()]

        class Tree:
            root = Root()

        tree = Tree()

        assert RobinsonFouldsDistance._get_all_bipartitions_direct(tree) is None
        assert RobinsonFouldsDistance._get_all_bipartition_id_sets_direct(tree) is None

    def test_calculate_robinson_foulds_distance_uses_cached_bipartitions(
        self, mocker, tree_simple, tree_simple_other, args
    ):
        rf = RobinsonFouldsDistance(args)
        mocker.patch.object(
            rf,
            "get_all_bipartitions",
            side_effect=AssertionError("standard trees should use compact id splits"),
        )
        mocker.patch.object(
            rf,
            "compare_trees_optimized",
            side_effect=AssertionError("fallback should not be used"),
        )

        plain_rf, normalized_rf = rf.calculate_robinson_foulds_distance(
            tree_simple, tree_simple_other
        )

        assert plain_rf == 8
        assert normalized_rf == pytest.approx(0.8)

    def test_calculate_robinson_foulds_distance_uses_direct_terminal_count(
        self, monkeypatch, args
    ):
        rf = RobinsonFouldsDistance(args)
        tree_zero = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        tree_one = copy.deepcopy(tree_zero)

        def fail_count_terminals(*_args, **_kwargs):
            raise AssertionError("generic terminal count should not be used")

        monkeypatch.setattr(type(tree_zero), "count_terminals", fail_count_terminals)

        plain_rf, normalized_rf = rf.calculate_robinson_foulds_distance(
            tree_zero,
            tree_one,
        )

        assert plain_rf == 0
        assert normalized_rf == 0

    def test_terminal_count_direct_counts_standard_tree_tips(self):
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);"), "newick")

        assert RobinsonFouldsDistance._terminal_count_direct(tree) == 5

    def test_first_terminal_name_uses_direct_tree_traversal(self, monkeypatch):
        tree = Phylo.read(StringIO("(((A:1,B:1):1,C:1):1,(D:1,E:1):1);"), "newick")

        def fail_get_terminals(self):
            raise AssertionError("get_terminals should not be used")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        assert RobinsonFouldsDistance._first_terminal_name(tree) == "A"

    def test_first_terminal_name_falls_back_for_nonstandard_clade_container(self):
        class NonstandardClade:
            name = "ignored"
            clades = ("not", "a", "list")

        class NonstandardTree:
            root = NonstandardClade()

            def get_terminals(self):
                terminal = type("Terminal", (), {"name": "fallback"})()
                return [terminal]

        assert RobinsonFouldsDistance._first_terminal_name(NonstandardTree()) == "fallback"

    def test_first_terminal_name_falls_back_for_inner_nonstandard_clade_container(self):
        class InnerNonstandardClade:
            name = "ignored"
            clades = ("not", "a", "list")

        class RootClade:
            clades = [InnerNonstandardClade()]

        class NonstandardTree:
            root = RootClade()

            def get_terminals(self):
                terminal = type("Terminal", (), {"name": "fallback"})()
                return [terminal]

        assert RobinsonFouldsDistance._first_terminal_name(NonstandardTree()) == "fallback"

    def test_calculate_multiple_rf_distances_sequential_path(self, mocker, args):
        rf = RobinsonFouldsDistance(args)

        tree_pairs = [("tree_zero", "tree_one") for _ in range(3)]
        expected = [(idx, idx / 10.0) for idx in range(len(tree_pairs))]
        mock_calc = mocker.patch.object(
            rf,
            "calculate_robinson_foulds_distance",
            side_effect=expected,
        )

        results = rf.calculate_multiple_rf_distances(tree_pairs)

        assert mock_calc.call_count == len(tree_pairs)
        assert results == expected

    def test_calculate_multiple_rf_distances_medium_input_skips_executor(self, mocker, args):
        rf = RobinsonFouldsDistance(args)

        tree_pairs = [("tree_zero", "tree_one") for _ in range(20)]
        expected = [(idx, idx / 10.0) for idx in range(len(tree_pairs))]
        mock_calc = mocker.patch.object(
            rf,
            "calculate_robinson_foulds_distance",
            side_effect=expected,
        )
        mocked_executor = mocker.patch("phykit.services.tree.rf_distance.ProcessPoolExecutor")

        results = rf.calculate_multiple_rf_distances(tree_pairs)

        assert mock_calc.call_count == len(tree_pairs)
        assert results == expected
        mocked_executor.assert_not_called()

    def test_calculate_multiple_rf_distances_parallel_path(self, mocker, args):
        rf = RobinsonFouldsDistance(args)
        rf.MP_MIN_TREE_PAIRS = 5

        class _DummyTree:
            def __init__(self, bipartitions, tip_count):
                self.bipartitions = bipartitions
                self._tip_count = tip_count

            def count_terminals(self):
                return self._tip_count

        created_executors = []

        class DummyExecutor:
            def __init__(self, max_workers):
                self.max_workers = max_workers
                created_executors.append(self)

            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, tb):
                return False

            def submit(self, fn, *args, **kwargs):
                future = Future()
                future.set_result(fn(*args, **kwargs))
                return future

        mocker.patch("phykit.services.tree.rf_distance.ProcessPoolExecutor", DummyExecutor)
        mocker.patch("phykit.services.tree.rf_distance.pickle.dumps", side_effect=lambda obj: obj)
        mocker.patch("phykit.services.tree.rf_distance.pickle.loads", side_effect=lambda obj: obj)
        mocker.patch.object(
            RobinsonFouldsDistance,
            "get_all_bipartitions",
            new=lambda self, tree: tree.bipartitions,
            create=True,
        )

        tree_pairs = [
            (
                _DummyTree({frozenset({1, 2}), frozenset({3, 4})}, 6),
                _DummyTree({frozenset({1, 2}), frozenset({4, 5})}, 6),
            ),
            (
                _DummyTree({frozenset({1, 2}), frozenset({2, 3})}, 6),
                _DummyTree({frozenset({1, 2}), frozenset({2, 3})}, 6),
            ),
            (
                _DummyTree({frozenset({1, 3})}, 5),
                _DummyTree({frozenset({2, 4})}, 5),
            ),
            (
                _DummyTree({frozenset({1, 4}), frozenset({2, 5})}, 7),
                _DummyTree({frozenset({1, 4})}, 7),
            ),
            (
                _DummyTree(set(), 8),
                _DummyTree({frozenset({1, 2})}, 8),
            ),
        ]

        results = rf.calculate_multiple_rf_distances(tree_pairs)

        assert created_executors

        expected = []
        for zero_tree, one_tree in tree_pairs:
            plain_rf = len(zero_tree.bipartitions ^ one_tree.bipartitions)
            normalized_rf = plain_rf / (2 * (zero_tree.count_terminals() - 3))
            expected.append((plain_rf, normalized_rf))

        for observed, anticipated in zip(results, expected):
            assert observed[0] == anticipated[0]
            assert observed[1] == pytest.approx(anticipated[1])

    def test_run_text_output_with_pruning(self, mocker, args):
        rf = RobinsonFouldsDistance(args)
        tree_zero = mocker.Mock()
        tree_one = mocker.Mock()
        term_zero = mocker.Mock()
        term_zero.name = "A"
        term_one = mocker.Mock()
        term_one.name = "A"
        tree_zero.get_terminals.return_value = [term_zero]
        tree_one.get_terminals.return_value = [term_one]
        mocker.patch.object(rf, "read_tree_file", return_value=tree_zero)
        mocker.patch.object(rf, "read_tree1_file", return_value=tree_one)
        mocker.patch.object(
            rf,
            "get_tip_names_from_tree",
            side_effect=[["A", "B"], ["A", "C"]],
        )
        mocker.patch.object(rf, "prune_tree_using_taxa_list", side_effect=[tree_zero, tree_one])
        mocker.patch.object(rf, "calculate_robinson_foulds_distance", return_value=(4, 0.5))
        mocked_print = mocker.patch("builtins.print")

        rf.run()

        assert rf.prune_tree_using_taxa_list.call_count == 2
        tree_zero.root_with_outgroup.assert_called_once_with("A")
        tree_one.root_with_outgroup.assert_called_once_with("A")
        mocked_print.assert_called_once_with("4\t0.5")

    def test_run_same_path_shortcuts_identical_rf(self, mocker):
        args = Namespace(tree_zero="/some/path/to/file.tre", tree_one="/some/path/to/file.tre")
        rf = RobinsonFouldsDistance(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        mocker.patch.object(rf, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(
            rf,
            "read_tree_file",
            side_effect=AssertionError("same-path RF should not copy tree zero"),
        )
        mocker.patch.object(
            rf,
            "read_tree1_file",
            side_effect=AssertionError("same-path RF should not copy tree one"),
        )
        mocker.patch.object(
            rf,
            "calculate_robinson_foulds_distance",
            side_effect=AssertionError("same-path RF should skip split comparison"),
        )
        mocker.patch.object(
            tree,
            "root_with_outgroup",
            side_effect=AssertionError("same-path RF should not mutate tree rooting"),
        )
        mocked_print = mocker.patch("builtins.print")

        rf.run()

        rf.read_tree_file_unmodified.assert_called_once_with()
        mocked_print.assert_called_once_with("0\t0.0")

    def test_run_json_output(self, mocker):
        args = Namespace(
            tree_zero="/some/path/to/file_zero.tre",
            tree_one="/some/path/to/file_one.tre",
            json=True,
        )
        rf = RobinsonFouldsDistance(args)
        tree_zero = mocker.Mock()
        tree_one = mocker.Mock()
        term_zero = mocker.Mock()
        term_zero.name = "A"
        term_one = mocker.Mock()
        term_one.name = "A"
        tree_zero.get_terminals.return_value = [term_zero]
        tree_one.get_terminals.return_value = [term_one]
        mocker.patch.object(rf, "read_tree_file", return_value=tree_zero)
        mocker.patch.object(rf, "read_tree1_file", return_value=tree_one)
        mocker.patch.object(rf, "get_tip_names_from_tree", side_effect=[["A"], ["A"]])
        mocker.patch.object(rf, "calculate_robinson_foulds_distance", return_value=(3, 0.375))
        mocked_json = mocker.patch("phykit.services.tree.rf_distance.print_json")

        rf.run()

        payload = mocked_json.call_args.args[0]
        assert payload["plain_rf"] == 3
        assert payload["normalized_rf"] == 0.375
