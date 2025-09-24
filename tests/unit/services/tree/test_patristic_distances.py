import pytest
from argparse import Namespace
from itertools import combinations
from math import isclose

from phykit.services.tree.patristic_distances import PatristicDistances


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

    def test_calculate_distance_between_pairs_parallel_path(self, mocker, args):
        t = PatristicDistances(args)
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
