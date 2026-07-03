import pytest
import subprocess
import sys
from argparse import Namespace
from concurrent.futures import Future
from itertools import combinations
from math import isclose
from Bio.Phylo.Newick import Clade, Tree

import phykit.services.tree.lb_score as lb_score_module
from phykit.services.tree.lb_score import LBScore


def test_module_import_defers_optional_progress_import():
    code = (
        "import sys; "
        "import phykit.services.tree.lb_score; "
        "assert 'typing' not in sys.modules; "
        "assert 'tqdm' not in sys.modules; "
        "assert 'Bio.Phylo' not in sys.modules; "
        "assert 'numpy' not in sys.modules; "
        "assert 'pickle' not in sys.modules; "
        "assert 'json' not in sys.modules; "
        "assert 'statistics' not in sys.modules; "
        "assert 'phykit.helpers.stats_summary' not in sys.modules; "
        "assert 'phykit.helpers.json_output' not in sys.modules; "
        "assert 'concurrent.futures' not in sys.modules; "
        "assert 'multiprocessing' not in sys.modules"
    )

    subprocess.run([sys.executable, "-c", code], check=True)


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

    def test_batched_tip_pairs_streams_all_pairs(self, args):
        t = LBScore(args)
        tips = ["tip0", "tip1", "tip2", "tip3", "tip4"]

        batches = list(t._batched_tip_pairs(tips, 3))
        observed_pairs = [pair for batch in batches for pair in batch]

        assert [len(batch) for batch in batches] == [3, 3, 3, 1]
        assert observed_pairs == list(combinations(tips, 2))

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

    def test_calculate_average_distance_between_tips_fast_path(self, mocker, args):
        t = LBScore(args)
        tips = ["tip0", "tip1", "tip2", "tip3"]
        tree = Tree(
            root=Clade(
                clades=[
                    Clade(branch_length=1.0, clades=[
                        Clade(branch_length=0.1, name="tip0"),
                        Clade(branch_length=0.2, name="tip1"),
                    ]),
                    Clade(branch_length=2.0, clades=[
                        Clade(branch_length=0.3, name="tip2"),
                        Clade(branch_length=0.4, name="tip3"),
                    ]),
                ],
            )
        )
        mocked_executor = mocker.patch("phykit.services.tree.lb_score.ProcessPoolExecutor")
        pairs = list(combinations(tips, 2))
        expected = sum(tree.distance(tip1, tip2) for tip1, tip2 in pairs) / len(pairs)

        result = t.calculate_average_distance_between_tips(tips, tree)

        assert result == pytest.approx(expected)
        mocked_executor.assert_not_called()

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

    def test_average_distance_to_other_taxa_reuses_tip_set(self, mocker, args):
        class CountingTips(list):
            def __init__(self, values):
                super().__init__(values)
                self.iterations = 0

            def __iter__(self):
                self.iterations += 1
                return super().__iter__()

        t = LBScore(args)
        mocker.patch.object(t, "calculate_pairwise_tip_distances_fast", return_value=None)
        tips = CountingTips([f"tip{i}" for i in range(12)])
        tree = _IndexedDummyTree()

        t.calculate_average_distance_of_taxon_to_other_taxa(tips, tree)

        assert tips.iterations == 2

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

    def test_calculate_average_distance_of_taxon_to_other_taxa_fast_path_preserves_bug(
        self, mocker, args
    ):
        t = LBScore(args)
        tips = ["tip0", "tip1", "tip2", "tip3"]
        tree = Tree(
            root=Clade(
                clades=[
                    Clade(branch_length=1.0, clades=[
                        Clade(branch_length=0.1, name="tip0"),
                        Clade(branch_length=0.2, name="tip1"),
                    ]),
                    Clade(branch_length=2.0, clades=[
                        Clade(branch_length=0.3, name="tip2"),
                        Clade(branch_length=0.4, name="tip3"),
                    ]),
                ],
            )
        )
        mocked_executor = mocker.patch("phykit.services.tree.lb_score.ProcessPoolExecutor")
        expected = []
        for tip in tips:
            tips_minus_i = list(set(tips) - set(tip))
            distances = [tree.distance(tip, other) for other in tips_minus_i]
            expected.append(sum(distances) / len(distances) if distances else 0)

        result = t.calculate_average_distance_of_taxon_to_other_taxa(tips, tree)

        assert result == pytest.approx(expected)
        mocked_executor.assert_not_called()

    def test_calculate_average_distance_of_taxon_to_other_taxa_fast_path_accumulates_pairwise_distances(
        self, args
    ):
        t = LBScore(args)
        tips = ["tip0", "tip1", "tip2"]
        pairwise_distances = (
            [("tip0", "tip1"), ("tip0", "tip2"), ("tip1", "tip2")],
            [3.0, 6.0, 9.0],
        )

        result = t.calculate_average_distance_of_taxon_to_other_taxa(
            tips,
            _IndexedDummyTree(),
            pairwise_distances,
        )

        assert result == pytest.approx([3.0, 4.0, 5.0])

    def test_calculate_lb_score_uses_linear_tree_components(self, mocker, args):
        t = LBScore(args)
        tips = ["tip0", "tip1", "tip2", "tip3"]
        tree = Tree(
            root=Clade(
                clades=[
                    Clade(branch_length=1.0, clades=[
                        Clade(branch_length=0.1, name="tip0"),
                        Clade(branch_length=0.2, name="tip1"),
                    ]),
                    Clade(branch_length=2.0, clades=[
                        Clade(branch_length=0.3, name="tip2"),
                        Clade(branch_length=0.4, name="tip3"),
                    ]),
                ],
            )
        )
        pairwise_spy = mocker.patch.object(
            t,
            "calculate_pairwise_tip_distances_fast",
            wraps=t.calculate_pairwise_tip_distances_fast,
        )

        result_tips, lb_scores = t.calculate_lb_score(tree)

        assert result_tips == tips
        assert len(lb_scores) == len(tips)
        assert pairwise_spy.call_count == 0

    def test_calculate_lb_components_fast_matches_pairwise_cache(self, args):
        t = LBScore(args)
        tips = ["tip0", "tip1", "tip2", "tip3"]
        tree = Tree(
            root=Clade(
                clades=[
                    Clade(branch_length=1.0, clades=[
                        Clade(branch_length=0.1, name="tip0"),
                        Clade(branch_length=0.2, name="tip1"),
                    ]),
                    Clade(branch_length=2.0, clades=[
                        Clade(branch_length=0.3, name="tip2"),
                        Clade(branch_length=0.4, name="tip3"),
                    ]),
                ],
            )
        )
        pairwise_distances = t.calculate_pairwise_tip_distances_fast(tree, tips)
        expected_avg_dist = t.calculate_average_distance_between_tips(
            tips,
            tree,
            pairwise_distances,
        )
        expected_avg_pdis = t.calculate_average_distance_of_taxon_to_other_taxa(
            tips,
            tree,
            pairwise_distances,
        )

        result = LBScore._calculate_lb_components_fast(tree, tips)

        assert result is not None
        avg_dist, avg_pdis = result
        assert avg_dist == pytest.approx(expected_avg_dist)
        assert avg_pdis == pytest.approx(expected_avg_pdis)

    def test_calculate_lb_components_fast_rejects_duplicate_tips(self):
        tree = Tree(
            root=Clade(
                clades=[
                    Clade(branch_length=1.0, name="tip0"),
                    Clade(branch_length=1.0, name="tip1"),
                ],
            )
        )

        assert LBScore._calculate_lb_components_fast(
            tree,
            ["tip0", "tip0"],
        ) is None

    def test_calculate_lb_components_fast_handles_mixed_child_counts_without_reversed(
        self, args
    ):
        class NoReversedList(list):
            def __reversed__(self):
                raise AssertionError("LB component traversal should push children explicitly")

        class DirectClade:
            def __init__(self, branch_length=None, name=None, clades=None):
                self.branch_length = branch_length
                self.name = name
                self.clades = NoReversedList(clades or [])

        tips = ["A", "B", "C", "D", "E", "F"]
        t = LBScore(args)
        reference_tree = Tree(
            root=Clade(
                clades=[
                    Clade(branch_length=1.0, name="A"),
                    Clade(
                        branch_length=2.0,
                        clades=[
                            Clade(branch_length=3.0, name="B"),
                            Clade(branch_length=4.0, name="C"),
                        ],
                    ),
                    Clade(
                        branch_length=5.0,
                        clades=[
                            Clade(branch_length=6.0, name="D"),
                            Clade(branch_length=7.0, name="E"),
                            Clade(branch_length=8.0, name="F"),
                        ],
                    ),
                ],
            )
        )
        direct_tree = type(
            "Tree",
            (),
            {
                "root": DirectClade(
                    clades=[
                        DirectClade(1.0, "A"),
                        DirectClade(
                            2.0,
                            clades=[
                                DirectClade(3.0, "B"),
                                DirectClade(4.0, "C"),
                            ],
                        ),
                        DirectClade(
                            5.0,
                            clades=[
                                DirectClade(6.0, "D"),
                                DirectClade(7.0, "E"),
                                DirectClade(8.0, "F"),
                            ],
                        ),
                    ],
                )
            },
        )()
        pairwise_distances = t.calculate_pairwise_tip_distances_fast(
            reference_tree,
            tips,
        )
        expected_avg_dist = t.calculate_average_distance_between_tips(
            tips,
            reference_tree,
            pairwise_distances,
        )
        expected_avg_pdis = t.calculate_average_distance_of_taxon_to_other_taxa(
            tips,
            reference_tree,
            pairwise_distances,
        )

        result = LBScore._calculate_lb_components_fast(direct_tree, tips)

        assert result is not None
        avg_dist, avg_pdis = result
        assert avg_dist == pytest.approx(expected_avg_dist)
        assert avg_pdis == pytest.approx(expected_avg_pdis)

    def test_historical_denominator_quirk_preserved_without_set_copy(self, args):
        t = LBScore(args)
        tips = ["A", "B", "AB", "BA"]
        tree = Tree(
            root=Clade(
                clades=[
                    Clade(branch_length=1.0, clades=[
                        Clade(branch_length=0.1, name="A"),
                        Clade(branch_length=0.2, name="B"),
                    ]),
                    Clade(branch_length=2.0, clades=[
                        Clade(branch_length=0.3, name="AB"),
                        Clade(branch_length=0.4, name="BA"),
                    ]),
                ],
            )
        )
        pairwise_distances = t.calculate_pairwise_tip_distances_fast(tree, tips)
        combos, distances = pairwise_distances
        distance_sums = {tip: 0.0 for tip in tips}
        for (tip_a, tip_b), distance in zip(combos, distances):
            distance_sums[tip_a] += distance
            distance_sums[tip_b] += distance
        tip_set = set(tips)
        expected_avg_pdis = [
            distance_sums[tip] / len(tip_set - set(tip))
            if len(tip_set - set(tip))
            else 0
            for tip in tips
        ]

        observed_pairwise = t.calculate_average_distance_of_taxon_to_other_taxa(
            tips,
            tree,
            pairwise_distances,
        )
        observed_fast = LBScore._calculate_lb_components_fast(tree, tips)[1]

        assert observed_pairwise == pytest.approx(expected_avg_pdis)
        assert observed_fast == pytest.approx(expected_avg_pdis)

    def test_historical_denominator_scan_preserves_duplicate_character_quirk(self):
        tip_count = lb_score_module._DENOMINATOR_SCAN_MIN_TIPS
        tip_set = {"A", "B", "AABB"}

        assert (
            LBScore._historical_other_taxa_denominator("AABB", tip_set, tip_count)
            == tip_count - 2
        )

    def test_historical_denominator_large_tip_count_avoids_set_constructor(
        self,
        monkeypatch,
    ):
        tip_count = lb_score_module._DENOMINATOR_SCAN_MIN_TIPS
        tip_set = {f"taxon_{idx}" for idx in range(tip_count)}

        def fail_set(*_args, **_kwargs):
            raise AssertionError("large denominator path should scan tip characters")

        monkeypatch.setattr(lb_score_module, "set", fail_set, raising=False)

        assert (
            LBScore._historical_other_taxa_denominator(
                "taxon_5000",
                tip_set,
                tip_count,
            )
            == tip_count
        )

    def test_calculate_lb_score_falls_back_for_nonstandard_tree(self, mocker, args):
        t = LBScore(args)
        tree = _IndexedDummyTree()
        tips = ["tip0", "tip1", "tip2"]
        mocker.patch.object(t, "get_tip_names_from_tree", return_value=tips)
        pairwise_spy = mocker.patch.object(
            t,
            "calculate_pairwise_tip_distances_fast",
            return_value=(
                [("tip0", "tip1"), ("tip0", "tip2"), ("tip1", "tip2")],
                [3.0, 6.0, 9.0],
            ),
        )

        result_tips, lb_scores = t.calculate_lb_score(tree)

        assert result_tips == tips
        assert len(lb_scores) == len(tips)
        pairwise_spy.assert_called_once_with(tree, tips)

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

    def test_calculate_lb_score_per_taxa_handles_fractional_values(self, args):
        t = LBScore(args)

        avg_pdis = [0.25, 0.5, 1.0]
        avg_dist = 0.5

        result = t.calculate_lb_score_per_taxa(avg_pdis, avg_dist)

        assert result == pytest.approx([-50.0, 0.0, 100.0])

    def test_calculate_lb_score_per_taxa_does_not_import_numpy(self):
        code = """
from argparse import Namespace
import sys
from phykit.services.tree.lb_score import LBScore
svc = LBScore(Namespace(tree="unused", verbose=False, json=False))
assert svc.calculate_lb_score_per_taxa([10.0, 20.0], 10.0) == [0.0, 100.0]
assert "numpy" not in sys.modules
"""
        subprocess.run([sys.executable, "-c", code], check=True)

    def test_run_verbose_batches_lb_score_rows(self, mocker, args, capsys):
        args.verbose = True
        args.json = False
        t = LBScore(args)
        mocker.patch.object(t, "read_tree_file_unmodified", return_value=object())
        mocker.patch.object(
            t,
            "calculate_lb_score",
            return_value=(["a", "b"], [1.23456, -2.0]),
        )

        t.run()

        captured = capsys.readouterr()
        assert captured.out == "a\t1.2346\nb\t-2.0\n"

    def test_run_verbose_empty_rows_prints_nothing(self, mocker, args, capsys):
        args.verbose = True
        args.json = False
        t = LBScore(args)
        mocker.patch.object(t, "read_tree_file_unmodified", return_value=object())
        mocker.patch.object(t, "calculate_lb_score", return_value=([], []))

        t.run()

        captured = capsys.readouterr()
        assert captured.out == ""

    def test_run_uses_unmodified_tree_read(self, mocker, args):
        args.verbose = False
        args.json = False
        tree = object()
        t = LBScore(args)
        read_tree = mocker.patch.object(
            t, "read_tree_file_unmodified", return_value=tree
        )
        mocker.patch.object(
            t,
            "read_tree_file",
            side_effect=AssertionError("copying tree reader should not be used"),
        )
        calc = mocker.patch.object(
            t,
            "calculate_lb_score",
            return_value=(["a", "b"], [1.0, 2.0]),
        )
        mocker.patch(
            "phykit.services.tree.lb_score.calculate_summary_statistics_from_arr",
            return_value={"mean": 1.5},
        )
        mocker.patch("phykit.services.tree.lb_score.print_summary_statistics")

        t.run()

        read_tree.assert_called_once_with()
        calc.assert_called_once_with(tree)
