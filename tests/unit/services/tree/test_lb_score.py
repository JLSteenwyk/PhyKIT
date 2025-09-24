import pytest
from argparse import Namespace
from concurrent.futures import Future
from itertools import combinations
from math import isclose

from phykit.services.tree.lb_score import LBScore


class _IndexedDummyTree:
    def _index(self, tip: str) -> int:
        digits = "".join(ch for ch in tip if ch.isdigit())
        return int(digits) if digits else 0

    def distance(self, tip1: str, tip2: str) -> int:
        return abs(self._index(tip1) - self._index(tip2))


@pytest.fixture
def args():
    kwargs = dict(tree="/some/path/to/file.tre", verbose=None)
    return Namespace(**kwargs)


class TestLBScore(object):
    def test_init_sets_tree_file_path(self, args):
        t = LBScore(args)
        assert t.tree_file_path == args.tree
        assert t.output_file_path is None

    def test_read_file_reads_tree_file_path(self, mocker, args):
        # Mock the cached tree read method instead of Phylo.read
        mock_cached_read = mocker.patch("phykit.services.tree.base.Tree._cached_tree_read")
        mock_get_hash = mocker.patch("phykit.services.tree.base.Tree._get_file_hash", return_value="test_hash")

        t = LBScore(args)
        t.read_tree_file()

        # Verify the cached read was called with the correct parameters
        mock_get_hash.assert_called_with(args.tree)
        mock_cached_read.assert_called_with(args.tree, "newick", "test_hash")

    # def test_calculate_lb_score_zero_branch_len(self, tree_zero_branch_length, args):
    #     t = LBScore(args)
    #     res = t.calculate_lb_score(tree_zero_branch_length)
    #     assert res is None

    def test_calculate_treeness(self, tree_simple, args):
        t = LBScore(args)
        tips, LBis = t.calculate_lb_score(tree_simple)
        expected_tips = [
            "raccoon",
            "bear",
            "sea_lion",
            "seal",
            "monkey",
            "cat",
            "weasel",
            "dog",
        ]
        expected_LBis = [
            -27.07902352846223,
            -39.283360704291205,
            -31.053612361805104,
            -31.04770664682598,
            65.67086344271493,
            12.796396820175548,
            -28.5329449372696,
            -21.470612084236517,
        ]
        assert tips == expected_tips
        for idx, value in enumerate(LBis):
            assert isclose(value, expected_LBis[idx], rel_tol=0.001)

    def test_calculate_average_distance_between_tips_single_tip(self, args):
        t = LBScore(args)
        tips = ["tip0"]
        tree = _IndexedDummyTree()

        assert t.calculate_average_distance_between_tips(tips, tree) == 0

    def test_calculate_average_distance_between_tips_small_dataset(self, args):
        t = LBScore(args)
        tips = ["tip0", "tip2", "tip5"]
        tree = _IndexedDummyTree()
        pairs = list(combinations(tips, 2))
        expected = sum(tree.distance(t1, t2) for t1, t2 in pairs) / len(pairs)

        result = t.calculate_average_distance_between_tips(tips, tree)

        assert result == pytest.approx(expected)

    def test_calculate_average_distance_between_tips_parallel_path(
        self, mocker, args
    ):
        t = LBScore(args)

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

        mocker.patch("phykit.services.tree.lb_score.mp.cpu_count", return_value=4)
        mocker.patch("phykit.services.tree.lb_score.pickle.dumps", side_effect=lambda obj: obj)
        mocker.patch("phykit.services.tree.lb_score.pickle.loads", side_effect=lambda obj: obj)
        mocker.patch(
            "phykit.services.tree.lb_score.ProcessPoolExecutor",
            DummyExecutor,
        )

        tree = _IndexedDummyTree()
        tips = [f"tip{i}" for i in range(15)]
        pairs = list(combinations(tips, 2))
        expected = sum(tree.distance(t1, t2) for t1, t2 in pairs) / len(pairs)

        result = t.calculate_average_distance_between_tips(tips, tree)

        assert created_executors
        assert result == pytest.approx(expected)

    def test_calculate_average_distance_of_taxon_to_other_taxa_small_dataset(
        self, args
    ):
        t = LBScore(args)

        tips = ["tip0", "tip12"]
        tree = _IndexedDummyTree()
        expected = []

        for tip in tips:
            other_tips = list(set(tips) - set(tip))
            distances = [tree.distance(tip, other) for other in other_tips]
            avg = sum(distances) / len(distances) if distances else 0
            expected.append(avg)

        result = t.calculate_average_distance_of_taxon_to_other_taxa(tips, tree)

        assert result == pytest.approx(expected)

    def test_calculate_average_distance_of_taxon_to_other_taxa_parallel_path(
        self, mocker, args
    ):
        t = LBScore(args)

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

        mocker.patch("phykit.services.tree.lb_score.mp.cpu_count", return_value=4)
        mocker.patch("phykit.services.tree.lb_score.pickle.dumps", side_effect=lambda obj: obj)
        mocker.patch("phykit.services.tree.lb_score.pickle.loads", side_effect=lambda obj: obj)
        mocker.patch(
            "phykit.services.tree.lb_score.ProcessPoolExecutor",
            DummyExecutor,
        )

        tips = [f"tip{i}" for i in range(60)]
        tree = _IndexedDummyTree()
        expected = []

        for tip in tips:
            other_tips = list(set(tips) - set(tip))
            distances = [tree.distance(tip, other) for other in other_tips]
            avg = sum(distances) / len(distances) if distances else 0
            expected.append(avg)

        result = t.calculate_average_distance_of_taxon_to_other_taxa(tips, tree)

        assert created_executors
        assert result == pytest.approx(expected)

    def test_calculate_lb_score_per_taxa_zero_avg_dist_exits(self, args):
        t = LBScore(args)

        with pytest.raises(SystemExit):
            t.calculate_lb_score_per_taxa([1.0], 0)

    def test_calculate_lb_score_per_taxa_returns_expected_values(self, args):
        t = LBScore(args)

        avg_pdis = [10.0, 20.0]
        avg_dist = 10.0

        result = t.calculate_lb_score_per_taxa(avg_pdis, avg_dist)

        assert result == pytest.approx([0.0, 100.0])
