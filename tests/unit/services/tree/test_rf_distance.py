import copy
import pytest
from argparse import Namespace
from concurrent.futures import Future

from phykit.services.tree.rf_distance import RobinsonFouldsDistance


@pytest.fixture
def args():
    kwargs = dict(
        tree_zero="/some/path/to/file.tre", tree_one="/some/path/to/file.tre",
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

    def test_calculate_multiple_rf_distances_parallel_path(self, mocker, args):
        rf = RobinsonFouldsDistance(args)

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
