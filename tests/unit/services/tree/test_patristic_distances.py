import pytest
import numpy as np
import subprocess
import sys
from argparse import Namespace
from itertools import combinations
from math import isclose
from Bio.Phylo.Newick import Clade, Tree

import phykit.services.tree.patristic_distances as patristic_distances_module
from phykit.services.tree.patristic_distances import PatristicDistances


def test_module_import_does_not_import_multiprocessing_or_pickle():
    code = """
import sys
import phykit.services.tree.patristic_distances as module
assert hasattr(module.mp, "Pool")
assert hasattr(module.pickle, "dumps")
assert "typing" not in sys.modules
assert "multiprocessing" not in sys.modules
assert "pickle" not in sys.modules
assert "json" not in sys.modules
assert "statistics" not in sys.modules
assert "phykit.helpers.stats_summary" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "Bio.Phylo" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_lazy_pickle_caches_resolved_attributes():
    lazy_pickle = patristic_distances_module._LazyPickle()

    blob = lazy_pickle.dumps({"value": 1})

    assert lazy_pickle._module is not None
    assert lazy_pickle.__dict__["dumps"] is lazy_pickle.dumps
    assert lazy_pickle.__dict__["loads"] is lazy_pickle.loads
    assert lazy_pickle.loads(blob) == {"value": 1}


def test_lazy_multiprocessing_caches_module_and_keeps_cpu_count_patchable(monkeypatch):
    import multiprocessing

    lazy_mp = patristic_distances_module._LazyMultiprocessing()

    monkeypatch.setattr(multiprocessing, "cpu_count", lambda: 11)
    assert lazy_mp.cpu_count() == 11

    monkeypatch.setattr(multiprocessing, "cpu_count", lambda: 7)
    assert lazy_mp.cpu_count() == 7
    assert lazy_mp._module is multiprocessing


@pytest.fixture
def args():
    kwargs = dict(tree="/some/path/to/file.tre", verbose=None)
    return Namespace(**kwargs)


class _IndexedDummyTree:
    def _index(self, tip: str) -> int:
        digits = "".join(ch for ch in tip if ch.isdigit())
        return int(digits) if digits else 0

    def distance(self, tip1: str, tip2: str) -> int:
        return abs(self._index(tip1) - self._index(tip2))


class TestPatristicDistances(object):
    def test_init_sets_tree_file_path(self, args):
        t = PatristicDistances(args)
        assert t.tree_file_path == args.tree
        assert t.output_file_path is None

    def test_read_file_reads_tree_file_path(self, mocker, args):
        # Mock the cached tree read method instead of Phylo.read
        mock_cached_read = mocker.patch("phykit.services.tree.base.Tree._cached_tree_read")
        mock_get_hash = mocker.patch("phykit.services.tree.base.Tree._get_file_hash", return_value="test_hash")

        t = PatristicDistances(args)
        t.read_tree_file()

        # Verify the cached read was called with the correct parameters
        mock_get_hash.assert_called_with(args.tree)
        mock_cached_read.assert_called_with(args.tree, "newick", "test_hash")

    def test_calculate_patristic_distances(self, tree_simple, args):
        t = PatristicDistances(args)
        patristic_distances, combos, stats = t.calculate_patristic_distances(
            tree_simple
        )
        assert isinstance(stats["mean"], float)
        assert isinstance(stats["median"], float)
        assert isinstance(stats["twenty_fifth"], float)
        assert isinstance(stats["seventy_fifth"], float)
        assert isinstance(stats["standard_deviation"], float)
        assert isinstance(stats["variance"], float)
        assert isinstance(stats["minimum"], float)
        assert isinstance(stats["maximum"], float)
        assert isclose(stats["mean"], 76.19737857142857, rel_tol=0.001)
        assert isclose(stats["median"], 49.588789999999996, rel_tol=0.001)
        assert isclose(stats["twenty_fifth"], 40.50536, rel_tol=0.001)
        assert isclose(stats["seventy_fifth"], 108.13853, rel_tol=0.001)
        assert isclose(stats["standard_deviation"], 45.46979239234539, rel_tol=0.001)
        assert isclose(stats["variance"], 2067.5020202029905, rel_tol=0.001)
        assert isclose(stats["minimum"], 24.0, rel_tol=0.001)
        assert isclose(stats["maximum"], 152.88127, rel_tol=0.001)

    def test_calculate_distance_between_pairs_small_dataset(self, args):
        t = PatristicDistances(args)
        tips = ["tip0", "tip2", "tip5"]
        tree = _IndexedDummyTree()

        combos, patristic_distances = t.calculate_distance_between_pairs(tips, tree)

        expected_combos = list(combinations(tips, 2))
        assert combos == expected_combos

        expected_distances = [tree.distance(*combo) for combo in expected_combos]
        assert patristic_distances == expected_distances

    def test_calculate_distance_batch_accepts_iterable_pairs(self, args):
        class Pair:
            def __init__(self, tip_a, tip_b):
                self.tip_a = tip_a
                self.tip_b = tip_b

            def __iter__(self):
                yield self.tip_a
                yield self.tip_b

            def __getitem__(self, _idx):
                raise AssertionError("distance batches should unpack pairs directly")

        t = PatristicDistances(args)
        tree_pickle = patristic_distances_module.pickle.dumps(_IndexedDummyTree())

        observed = t._calculate_distance_batch(
            tree_pickle,
            [Pair("tip0", "tip2"), Pair("tip2", "tip5")],
        )

        assert observed == [2, 3]

    def test_calculate_distance_between_pairs_fast_matches_biopython(self, mocker, args):
        t = PatristicDistances(args)
        tips = [f"tip{i}" for i in range(12)]
        tree = Tree(
            root=Clade(
                clades=[
                    Clade(branch_length=1.0, clades=[
                        Clade(branch_length=0.1, name="tip0"),
                        Clade(branch_length=0.2, name="tip1"),
                        Clade(branch_length=0.3, name="tip2"),
                        Clade(branch_length=0.4, name="tip3"),
                        Clade(branch_length=0.5, name="tip4"),
                        Clade(branch_length=0.6, name="tip5"),
                    ]),
                    Clade(branch_length=2.0, clades=[
                        Clade(branch_length=0.7, name="tip6"),
                        Clade(branch_length=0.8, name="tip7"),
                        Clade(branch_length=0.9, name="tip8"),
                        Clade(branch_length=1.0, name="tip9"),
                        Clade(branch_length=1.1, name="tip10"),
                        Clade(branch_length=1.2, name="tip11"),
                    ]),
                ],
            )
        )
        mocked_pool = mocker.patch("phykit.services.tree.patristic_distances.mp.Pool")

        combos, patristic_distances = t.calculate_distance_between_pairs(tips, tree)

        expected_combos = list(combinations(tips, 2))
        expected_distances = [tree.distance(*combo) for combo in expected_combos]
        assert combos == expected_combos
        assert np.allclose(patristic_distances, expected_distances)
        mocked_pool.assert_not_called()

    def test_calculate_distance_values_between_pairs_fast_skips_combos(self, mocker, args):
        t = PatristicDistances(args)
        tips = [f"tip{i}" for i in range(6)]
        tree = Tree(
            root=Clade(
                clades=[
                    Clade(branch_length=1.0, clades=[
                        Clade(branch_length=0.1, name="tip0"),
                        Clade(branch_length=0.2, name="tip1"),
                        Clade(branch_length=0.3, name="tip2"),
                    ]),
                    Clade(branch_length=2.0, clades=[
                        Clade(branch_length=0.4, name="tip3"),
                        Clade(branch_length=0.5, name="tip4"),
                        Clade(branch_length=0.6, name="tip5"),
                    ]),
                ],
            )
        )
        mocker.patch.object(
            PatristicDistances,
            "calculate_distance_between_pairs",
            side_effect=AssertionError("stats-only path should not build combos"),
        )

        patristic_distances = t.calculate_distance_values_between_pairs(tips, tree)

        expected_distances = [
            tree.distance(*combo)
            for combo in combinations(tips, 2)
        ]
        assert np.allclose(patristic_distances, expected_distances)

    def test_batched_tip_pairs_streams_all_pairs(self, args):
        t = PatristicDistances(args)
        tips = ["tip0", "tip1", "tip2", "tip3", "tip4"]

        batches = list(t._batched_tip_pairs(tips, 3))
        observed_pairs = [pair for batch in batches for pair in batch]

        assert [len(batch) for batch in batches] == [3, 3, 3, 1]
        assert observed_pairs == list(combinations(tips, 2))

    def test_distance_values_fallback_skips_returned_combos(self, mocker, args):
        t = PatristicDistances(args)
        tips = [f"tip{i}" for i in range(15)]
        tree = _IndexedDummyTree()

        class DummyPool:
            def __init__(self, processes):
                self.processes = processes

            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, tb):
                return False

            def map(self, func, chunks):
                return [func(chunk) for chunk in chunks]

            def imap(self, func, chunks):
                for chunk in chunks:
                    yield func(chunk)

        mocker.patch.object(
            PatristicDistances,
            "calculate_distance_between_pairs",
            side_effect=AssertionError("stats-only fallback should not return combos"),
        )
        mocker.patch("phykit.services.tree.patristic_distances.mp.Pool", DummyPool)
        mocker.patch("phykit.services.tree.patristic_distances.mp.cpu_count", return_value=4)
        mocker.patch("phykit.services.tree.patristic_distances.pickle.dumps", side_effect=lambda obj: obj)
        mocker.patch("phykit.services.tree.patristic_distances.pickle.loads", side_effect=lambda obj: obj)
        mocker.patch("phykit.services.tree.patristic_distances.sys.stderr.isatty", return_value=False)

        patristic_distances = t.calculate_distance_values_between_pairs(tips, tree)

        expected_distances = [
            tree.distance(*combo)
            for combo in combinations(tips, 2)
        ]
        assert patristic_distances == expected_distances

    def test_distance_values_fast_handles_mixed_child_counts_without_distance(
        self, monkeypatch, args
    ):
        from Bio.Phylo.BaseTree import TreeMixin

        t = PatristicDistances(args)
        tree = Tree(
            root=Clade(
                clades=[
                    Clade(
                        branch_length=1.0,
                        clades=[
                            Clade(branch_length=0.5, name="tip0"),
                            Clade(branch_length=0.7, name="tip1"),
                            Clade(branch_length=0.9, name="tip2"),
                        ],
                    ),
                    Clade(
                        branch_length=2.0,
                        clades=[
                            Clade(branch_length=1.1, name="tip3"),
                            Clade(
                                branch_length=1.3,
                                clades=[
                                    Clade(branch_length=0.2, name="tip4"),
                                    Clade(branch_length=0.4, name="tip5"),
                                ],
                            ),
                        ],
                    ),
                ],
            )
        )
        tips = [f"tip{i}" for i in range(6)]
        expected_distances = [
            tree.distance(*combo)
            for combo in combinations(tips, 2)
        ]

        def fail_distance(*_args, **_kwargs):
            raise AssertionError("fast path should use cached depths")

        monkeypatch.setattr(TreeMixin, "distance", fail_distance)

        patristic_distances = t.calculate_pairwise_tip_distance_values_fast(
            tree, tips
        )

        assert np.allclose(patristic_distances, expected_distances)

    def test_deep_tree_distance_values_use_lca_index_path(self, monkeypatch, args):
        t = PatristicDistances(args)
        root = Clade()
        current = root
        tips = []
        for idx in range(7):
            tip_name = f"tip{idx}"
            tips.append(tip_name)
            next_node = Clade(branch_length=1.0)
            current.clades = [
                Clade(branch_length=1.0, name=tip_name),
                next_node,
            ]
            current = next_node
        tips.append("tip7")
        current.clades = [Clade(branch_length=1.0, name="tip7")]
        tree = Tree(root=root)

        def fail_path_helper(*_args, **_kwargs):
            raise AssertionError("deep trees should avoid materialized tip paths")

        monkeypatch.setattr(PatristicDistances, "_PAIRWISE_LCA_DEPTH_THRESHOLD", 4)
        monkeypatch.setattr(
            PatristicDistances,
            "_pairwise_tip_distance_values_from_paths",
            staticmethod(fail_path_helper),
        )

        observed = t.calculate_pairwise_tip_distance_values_fast(tree, tips)
        expected = [
            tree.distance(*combo)
            for combo in combinations(tips, 2)
        ]

        assert np.allclose(observed, expected)

    def test_calculate_patristic_distance_stats_matches_full_stats(self, tree_simple, args):
        t = PatristicDistances(args)

        _, _, full_stats = t.calculate_patristic_distances(tree_simple)
        stats_only = t.calculate_patristic_distance_stats(tree_simple)

        assert stats_only == full_stats

    def test_calculate_distance_between_pairs_medium_fallback_skips_pool(self, mocker, args):
        t = PatristicDistances(args)
        tips = [f"tip{i}" for i in range(20)]
        tree = _IndexedDummyTree()
        mocked_pool = mocker.patch("phykit.services.tree.patristic_distances.mp.Pool")

        combos, patristic_distances = t.calculate_distance_between_pairs(tips, tree)

        expected_combos = list(combinations(tips, 2))
        expected_distances = [tree.distance(*combo) for combo in expected_combos]
        assert combos == expected_combos
        assert patristic_distances == expected_distances
        mocked_pool.assert_not_called()

    def test_distance_values_medium_fallback_skips_pool(self, mocker, args):
        t = PatristicDistances(args)
        tips = [f"tip{i}" for i in range(20)]
        tree = _IndexedDummyTree()
        mocked_pool = mocker.patch("phykit.services.tree.patristic_distances.mp.Pool")

        patristic_distances = t.calculate_distance_values_between_pairs(tips, tree)

        expected_distances = [
            tree.distance(*combo)
            for combo in combinations(tips, 2)
        ]
        assert patristic_distances == expected_distances
        mocked_pool.assert_not_called()

    def test_calculate_distance_between_pairs_parallel_path(self, mocker, args):
        t = PatristicDistances(args)
        t.MP_MIN_PAIRS = 100
        tips = [f"tip{i}" for i in range(15)]  # generates >100 combinations
        tree = _IndexedDummyTree()

        created_pools = []
        recorded_chunks = []

        class DummyPool:
            def __init__(self, processes):
                self.processes = processes
                created_pools.append(self)

            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, tb):
                return False

            def map(self, func, chunks):
                recorded_chunks.extend(chunks)
                return [func(chunk) for chunk in chunks]

            def imap(self, func, chunks):
                for chunk in chunks:
                    recorded_chunks.append(chunk)
                    yield func(chunk)

        mocker.patch("phykit.services.tree.patristic_distances.mp.Pool", DummyPool)
        mocker.patch("phykit.services.tree.patristic_distances.mp.cpu_count", return_value=4)
        mocker.patch("phykit.services.tree.patristic_distances.pickle.dumps", side_effect=lambda obj: obj)
        mocker.patch("phykit.services.tree.patristic_distances.pickle.loads", side_effect=lambda obj: obj)
        mocker.patch("phykit.services.tree.patristic_distances.sys.stderr.isatty", return_value=False)

        combos, patristic_distances = t.calculate_distance_between_pairs(tips, tree)

        assert created_pools
        assert all(len(chunk) > 0 for chunk in recorded_chunks)

        expected_combos = list(combinations(tips, 2))
        expected_distances = [tree.distance(*combo) for combo in expected_combos]

        assert combos == expected_combos
        assert patristic_distances == expected_distances

    def test_calculate_distance_between_pairs_parallel_tty_path(self, mocker, args):
        t = PatristicDistances(args)
        t.MP_MIN_PAIRS = 100
        tips = [f"tip{i}" for i in range(15)]
        tree = _IndexedDummyTree()

        class DummyPool:
            def __init__(self, processes):
                self.processes = processes

            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, tb):
                return False

            def map(self, func, chunks):
                return [func(chunk) for chunk in chunks]

            def imap(self, func, chunks):
                for chunk in chunks:
                    yield func(chunk)

        mocker.patch("phykit.services.tree.patristic_distances.mp.Pool", DummyPool)
        mocker.patch("phykit.services.tree.patristic_distances.mp.cpu_count", return_value=4)
        mocker.patch("phykit.services.tree.patristic_distances.pickle.dumps", side_effect=lambda obj: obj)
        mocker.patch("phykit.services.tree.patristic_distances.pickle.loads", side_effect=lambda obj: obj)
        mocker.patch("phykit.services.tree.patristic_distances.sys.stderr.isatty", return_value=True)
        mocked_progress = mocker.patch(
            "phykit.services.tree.patristic_distances._with_optional_progress",
            side_effect=lambda iterable, **kwargs: iterable,
        )

        combos, patristic_distances = t.calculate_distance_between_pairs(tips, tree)

        assert mocked_progress.called
        assert len(combos) == len(patristic_distances)

    def test_run_json_verbose(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", verbose=True, json=True)
        t = PatristicDistances(args)
        mocker.patch.object(t, "read_tree_file_unmodified", return_value=object())
        mocker.patch.object(
            t,
            "calculate_patristic_distances",
            return_value=([1.23456], [("a", "b")], {"mean": 1.23456}),
        )
        mocked_json = mocker.patch("phykit.services.tree.patristic_distances.print_json")

        t.run()

        payload = mocked_json.call_args.args[0]
        assert payload["verbose"] is True
        assert payload["rows"] == payload["pairs"]
        assert payload["rows"][0]["patristic_distance"] == 1.2346

    def test_run_json_verbose_uses_fast_rows(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", verbose=True, json=True)
        t = PatristicDistances(args)
        tree = Tree(
            root=Clade(
                clades=[
                    Clade(branch_length=1.0, name="a"),
                    Clade(branch_length=2.0, name="b"),
                    Clade(branch_length=3.0, name="c"),
                ],
            )
        )
        mocker.patch.object(t, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(
            t,
            "calculate_patristic_distances",
            side_effect=AssertionError("verbose JSON should use fast rows"),
        )
        mocked_json = mocker.patch("phykit.services.tree.patristic_distances.print_json")

        t.run()

        payload = mocked_json.call_args.args[0]
        assert payload["verbose"] is True
        assert payload["rows"] == [
            {"taxon_a": "a", "taxon_b": "b", "patristic_distance": 3.0},
            {"taxon_a": "a", "taxon_b": "c", "patristic_distance": 4.0},
            {"taxon_a": "b", "taxon_b": "c", "patristic_distance": 5.0},
        ]
        assert payload["pairs"] == payload["rows"]

    def test_run_json_non_verbose(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", verbose=False, json=True)
        t = PatristicDistances(args)
        mocker.patch.object(t, "read_tree_file_unmodified", return_value=object())
        mocker.patch.object(
            t,
            "calculate_patristic_distance_stats",
            return_value={"mean": 1.0},
        )
        mocked_json = mocker.patch("phykit.services.tree.patristic_distances.print_json")

        t.run()

        payload = mocked_json.call_args.args[0]
        assert payload["verbose"] is False
        assert payload["summary"]["mean"] == 1.0

    def test_run_verbose_handles_broken_pipe(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", verbose=True, json=False)
        t = PatristicDistances(args)
        mocker.patch.object(t, "read_tree_file_unmodified", return_value=object())
        mocker.patch.object(
            t,
            "calculate_patristic_distances",
            return_value=([1.0], [("a", "b")], {"mean": 1.0}),
        )
        mocker.patch("builtins.print", side_effect=BrokenPipeError())

        t.run()

    def test_run_verbose_batches_pairwise_rows(self, mocker, capsys):
        args = Namespace(tree="/some/path/to/file.tre", verbose=True, json=False)
        t = PatristicDistances(args)
        mocker.patch.object(t, "read_tree_file_unmodified", return_value=object())
        mocker.patch.object(
            t,
            "calculate_patristic_distances",
            return_value=(
                [1.23456, 2.0],
                [("a", "b"), ("a", "c")],
                {"mean": 1.61728},
            ),
        )

        t.run()

        captured = capsys.readouterr()
        assert captured.out == "a\tb\t1.2346\na\tc\t2.0\n"

    def test_run_verbose_text_uses_fast_rows(self, mocker, capsys):
        args = Namespace(tree="/some/path/to/file.tre", verbose=True, json=False)
        t = PatristicDistances(args)
        tree = Tree(
            root=Clade(
                clades=[
                    Clade(branch_length=1.0, name="a"),
                    Clade(branch_length=2.0, name="b"),
                    Clade(branch_length=3.0, name="c"),
                ],
            )
        )
        mocker.patch.object(t, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(
            t,
            "calculate_patristic_distances",
            side_effect=AssertionError("verbose text should use fast rows"),
        )

        t.run()

        captured = capsys.readouterr()
        assert captured.out == "a\tb\t3.0\na\tc\t4.0\nb\tc\t5.0\n"

    def test_run_verbose_empty_pairs_prints_nothing(self, mocker, capsys):
        args = Namespace(tree="/some/path/to/file.tre", verbose=True, json=False)
        t = PatristicDistances(args)
        mocker.patch.object(t, "read_tree_file_unmodified", return_value=object())
        mocker.patch.object(
            t,
            "calculate_patristic_distances",
            return_value=([], [], {}),
        )

        t.run()

        captured = capsys.readouterr()
        assert captured.out == ""

    def test_run_non_verbose_prints_summary(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", verbose=False, json=False)
        t = PatristicDistances(args)
        mocker.patch.object(t, "read_tree_file_unmodified", return_value=object())
        mocker.patch.object(
            t,
            "calculate_patristic_distance_stats",
            return_value={"mean": 1.0},
        )
        mocked_summary = mocker.patch("phykit.services.tree.patristic_distances.print_summary_statistics")

        t.run()

        mocked_summary.assert_called_once_with({"mean": 1.0})

    def test_run_uses_unmodified_tree_read(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", verbose=False, json=False)
        t = PatristicDistances(args)
        read_tree = mocker.patch.object(
            t, "read_tree_file_unmodified", return_value=object()
        )
        mocker.patch.object(
            t,
            "read_tree_file",
            side_effect=AssertionError("copying tree reader should not be used"),
        )
        mocker.patch.object(
            t,
            "calculate_patristic_distance_stats",
            return_value={"mean": 1.0},
        )
        mocker.patch("phykit.services.tree.patristic_distances.print_summary_statistics")

        t.run()

        read_tree.assert_called_once_with()
