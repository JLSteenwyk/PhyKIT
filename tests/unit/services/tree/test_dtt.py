"""
Unit tests for disparity through time (DTT).

Computes DTT curves and MDI following Harmon et al. (2003).
Cross-validated against R's geiger::dtt().

See tests/r_validation/validate_dtt.R for the R validation script.
"""
import pytest
import json
import os
import subprocess
import sys
import numpy as np
from argparse import Namespace
from io import StringIO
from pathlib import Path

from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin

from phykit.services.tree.dtt import Dtt, _max_terminal_depth, _mean_selected_columns
import phykit.services.tree.dtt as dtt_module


def test_mean_selected_columns_accumulates_small_selections(monkeypatch):
    values = np.arange(30, dtype=float).reshape(5, 6)

    def fail_mean(*_args, **_kwargs):
        raise AssertionError("small selected-column means should avoid np.mean")

    monkeypatch.setattr(dtt_module.np, "mean", fail_mean)

    np.testing.assert_allclose(
        _mean_selected_columns(values, [1]),
        values[:, 1],
    )
    np.testing.assert_allclose(
        _mean_selected_columns(values, [1, 3, 5]),
        values[:, [1, 3, 5]].mean(axis=1),
    )


def test_mean_selected_columns_uses_numpy_for_wide_selections():
    values = np.arange(120, dtype=float).reshape(10, 12)
    positions = list(range(10))

    np.testing.assert_allclose(
        _mean_selected_columns(values, positions),
        values[:, positions].mean(axis=1),
    )


def test_max_terminal_depth_scans_terminals_once():
    class SinglePassTerminals:
        def __init__(self, terminals):
            self.terminals = terminals
            self.iterations = 0

        def __iter__(self):
            self.iterations += 1
            if self.iterations > 1:
                raise AssertionError("terminals should be scanned once")
            return iter(self.terminals)

    terminals = [object(), object(), object()]
    depths = {terminals[0]: 1.5, terminals[1]: 4.0, terminals[2]: 2.25}
    iterable = SinglePassTerminals(terminals)

    assert _max_terminal_depth(iterable, depths, 0.5) == 3.5
    assert iterable.iterations == 1


def test_permutation_p_value_ge_counts_extreme_permutations(monkeypatch):
    permutations = np.array([0.5, 1.5, 2.5, 3.5])
    original_count_nonzero = dtt_module.np.count_nonzero
    calls = []

    def counting_count_nonzero(values):
        calls.append(values.copy())
        return original_count_nonzero(values)

    monkeypatch.setattr(dtt_module.np, "count_nonzero", counting_count_nonzero)

    p_value = dtt_module._permutation_p_value_ge(permutations, 2.0)

    assert p_value == pytest.approx(0.5)
    assert len(calls) == 1
    np.testing.assert_array_equal(calls[0], np.array([False, False, True, True]))
    assert np.isnan(dtt_module._permutation_p_value_ge(np.array([]), 2.0))


def test_batch_row_sum_squares_matches_explicit_squared_reduction():
    values = np.arange(2 * 3 * 4, dtype=float).reshape(2, 3, 4)

    observed = dtt_module._batch_row_sum_squares(values)
    expected = np.sum(values * values, axis=(1, 2))

    np.testing.assert_allclose(observed, expected)


def test_batch_trait_sums_use_array_method_for_multi_trait_cubes(monkeypatch):
    values = np.arange(2 * 3 * 4, dtype=float).reshape(2, 3, 4)

    def fail_sum(*_args, **_kwargs):
        raise AssertionError("multi-trait batch sums should use ndarray.sum")

    monkeypatch.setattr(dtt_module.np, "sum", fail_sum)

    np.testing.assert_allclose(
        dtt_module._batch_trait_sums(values),
        values.sum(axis=1),
    )


def test_batch_trait_sums_keep_numpy_path_for_single_trait_cubes(monkeypatch):
    values = np.arange(2 * 3, dtype=float).reshape(2, 3, 1)
    calls = []
    original_sum = dtt_module.np.sum

    def tracking_sum(*args, **kwargs):
        calls.append((args, kwargs))
        return original_sum(*args, **kwargs)

    monkeypatch.setattr(dtt_module.np, "sum", tracking_sum)

    np.testing.assert_allclose(
        dtt_module._batch_trait_sums(values),
        np.sum(values, axis=1),
    )
    assert len(calls) == 1


def test_batch_trait_sums_keep_numpy_path_for_wide_trait_cubes(monkeypatch):
    values = np.arange(2 * 3 * 8, dtype=float).reshape(2, 3, 8)
    calls = []
    original_sum = dtt_module.np.sum

    def tracking_sum(*args, **kwargs):
        calls.append((args, kwargs))
        return original_sum(*args, **kwargs)

    monkeypatch.setattr(dtt_module.np, "sum", tracking_sum)

    np.testing.assert_allclose(
        dtt_module._batch_trait_sums(values),
        np.sum(values, axis=1),
    )
    assert len(calls) == 1


def test_module_import_does_not_import_numpy():
    code = """
import sys
import phykit.services.tree.dtt as module
assert hasattr(module.np, "__getattr__")
assert callable(module.print_json)
assert callable(module.parse_multi_trait_file)
assert callable(module.trait_column_from_rows)
assert callable(module.trait_matrix_from_rows)
assert "json" not in sys.modules
assert "typing" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
assert "phykit.helpers.trait_parsing" not in sys.modules
assert "numpy" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_lazy_numpy_caches_resolved_attributes():
    lazy_np = dtt_module._LazyNumpy()

    count_nonzero_attr = lazy_np.count_nonzero

    assert lazy_np.__dict__["count_nonzero"] is count_nonzero_attr
    assert lazy_np.count_nonzero is count_nonzero_attr
    assert lazy_np._module is not None


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
ULTRAMETRIC_TREE = str(SAMPLE_FILES / "ultrametric_tree.tre")
MULTI_TRAITS = str(SAMPLE_FILES / "tree_simple_multi_traits.tsv")


def _make_args(**kwargs):
    defaults = dict(
        tree=ULTRAMETRIC_TREE,
        traits=MULTI_TRAITS,
        trait=None,
        index="avg_sq",
        nsim=0,
        seed=None,
        plot_output=None,
        json=False,
        fig_width=None,
        fig_height=None,
        dpi=300,
        no_title=False,
        title=None,
        legend_position=None,
        ylabel_fontsize=None,
        xlabel_fontsize=None,
        title_fontsize=None,
        axis_fontsize=None,
        colors=None,
        ladderize=False,
        cladogram=False,
        circular=False,
        color_file=None,
    )
    defaults.update(kwargs)
    return Namespace(**defaults)


class TestDttInit:
    def test_init_sets_fields(self):
        args = _make_args()
        svc = Dtt(args)
        assert svc.tree_file_path == ULTRAMETRIC_TREE
        assert svc.trait_data_path == MULTI_TRAITS
        assert svc.trait_column is None
        assert svc.index == "avg_sq"
        assert svc.nsim == 0
        assert svc.seed is None
        assert svc.plot_output is None
        assert svc.json_output is False

    def test_init_with_options(self):
        args = _make_args(trait="body_mass", index="avg_manhattan", nsim=100, seed=42)
        svc = Dtt(args)
        assert svc.trait_column == "body_mass"
        assert svc.index == "avg_manhattan"
        assert svc.nsim == 100
        assert svc.seed == 42


class TestComputeDisparity:
    def test_single_point_returns_zero(self):
        args = _make_args()
        svc = Dtt(args)
        data = np.array([[1.0, 2.0]])
        assert svc._compute_disparity(data) == 0.0

    def test_two_points_avg_sq(self):
        args = _make_args(index="avg_sq")
        svc = Dtt(args)
        data = np.array([[0.0], [3.0]])
        # (3-0)^2 = 9, one pair, avg = 9
        assert svc._compute_disparity(data) == pytest.approx(9.0)

    def test_two_points_avg_manhattan(self):
        args = _make_args(index="avg_manhattan")
        svc = Dtt(args)
        data = np.array([[0.0, 0.0], [3.0, 4.0]])
        # |3-0| + |4-0| = 7, one pair, avg = 7
        assert svc._compute_disparity(data) == pytest.approx(7.0)

    def test_three_points_avg_sq(self):
        args = _make_args(index="avg_sq")
        svc = Dtt(args)
        data = np.array([[0.0], [1.0], [3.0]])
        # pairs: (0,1)=1, (0,3)=9, (1,3)=4; avg = 14/3
        assert svc._compute_disparity(data) == pytest.approx(14.0 / 3.0)

    def test_avg_sq_disparity_uses_dot_sum_of_squares(self, monkeypatch):
        args = _make_args(index="avg_sq")
        svc = Dtt(args)
        data = np.array([[0.0, 2.0], [3.0, 4.0], [6.0, 8.0]])
        dot_calls = []
        original_dot = dtt_module.np.dot

        def counting_dot(left, right):
            dot_calls.append((left.shape, right.shape))
            return original_dot(left, right)

        monkeypatch.setattr(dtt_module.np, "dot", counting_dot)
        monkeypatch.setattr(
            dtt_module.np,
            "sum",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("avg_sq disparity should use dot reductions")
            ),
        )

        assert svc._compute_disparity(data) == pytest.approx(110.0 / 3.0)
        assert dot_calls == [((6,), (6,)), ((2,), (2,))]

    def test_three_points_avg_manhattan_multidimensional(self):
        args = _make_args(index="avg_manhattan")
        svc = Dtt(args)
        data = np.array([[0.0, 0.0], [1.0, 3.0], [4.0, 4.0]])
        # pairs: 4, 8, 4; avg = 16/3
        assert svc._compute_disparity(data) == pytest.approx(16.0 / 3.0)

    def test_avg_manhattan_disparity_uses_einsum_reduction(self, monkeypatch):
        args = _make_args(index="avg_manhattan")
        svc = Dtt(args)
        data = np.array([[0.0, 0.0], [1.0, 3.0], [4.0, 4.0]])

        monkeypatch.setattr(
            dtt_module.np,
            "sum",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("avg_manhattan disparity should use einsum")
            ),
        )

        assert svc._compute_disparity(data) == pytest.approx(16.0 / 3.0)

    def test_identical_values_return_zero(self):
        args = _make_args()
        svc = Dtt(args)
        data = np.array([[5.0], [5.0], [5.0]])
        assert svc._compute_disparity(data) == pytest.approx(0.0)


class TestComputeDtt:
    def test_dtt_starts_at_one(self):
        """DTT curve should start with relative disparity = 1.0 at time 0."""
        args = _make_args()
        svc = Dtt(args)
        tree = Phylo.read(ULTRAMETRIC_TREE, "newick")
        ordered_names = sorted(
            ["bear", "cat", "dog", "monkey", "raccoon", "sea_lion", "seal", "weasel"]
        )
        traits_data = {
            "raccoon": [1.04], "bear": [2.39], "sea_lion": [2.30], "seal": [1.88],
            "monkey": [0.60], "cat": [0.56], "weasel": [-0.30], "dog": [1.18],
        }
        data = np.array([traits_data[n] for n in ordered_names])
        times, dtt_values = svc._compute_dtt(tree, data, ordered_names)
        assert times[0] == pytest.approx(0.0)
        assert dtt_values[0] == pytest.approx(1.0)

    def test_dtt_ends_near_zero(self):
        """DTT curve should end with relative disparity near 0."""
        args = _make_args()
        svc = Dtt(args)
        tree = Phylo.read(ULTRAMETRIC_TREE, "newick")
        ordered_names = sorted(
            ["bear", "cat", "dog", "monkey", "raccoon", "sea_lion", "seal", "weasel"]
        )
        traits_data = {
            "raccoon": [1.04], "bear": [2.39], "sea_lion": [2.30], "seal": [1.88],
            "monkey": [0.60], "cat": [0.56], "weasel": [-0.30], "dog": [1.18],
        }
        data = np.array([traits_data[n] for n in ordered_names])
        times, dtt_values = svc._compute_dtt(tree, data, ordered_names)
        # The last DTT value should be 0 or very close to 0
        assert dtt_values[-1] == pytest.approx(0.0, abs=0.01)

    def test_times_monotonically_increasing(self):
        """Time points should be non-decreasing."""
        args = _make_args()
        svc = Dtt(args)
        tree = Phylo.read(ULTRAMETRIC_TREE, "newick")
        ordered_names = sorted(
            ["bear", "cat", "dog", "monkey", "raccoon", "sea_lion", "seal", "weasel"]
        )
        traits_data = {
            "raccoon": [1.04], "bear": [2.39], "sea_lion": [2.30], "seal": [1.88],
            "monkey": [0.60], "cat": [0.56], "weasel": [-0.30], "dog": [1.18],
        }
        data = np.array([traits_data[n] for n in ordered_names])
        times, _ = svc._compute_dtt(tree, data, ordered_names)
        for i in range(1, len(times)):
            assert times[i] >= times[i - 1]

    def test_dtt_values_non_negative(self):
        """Relative disparity should be >= 0."""
        args = _make_args()
        svc = Dtt(args)
        tree = Phylo.read(ULTRAMETRIC_TREE, "newick")
        ordered_names = sorted(
            ["bear", "cat", "dog", "monkey", "raccoon", "sea_lion", "seal", "weasel"]
        )
        traits_data = {
            "raccoon": [1.04], "bear": [2.39], "sea_lion": [2.30], "seal": [1.88],
            "monkey": [0.60], "cat": [0.56], "weasel": [-0.30], "dog": [1.18],
        }
        data = np.array([traits_data[n] for n in ordered_names])
        _, dtt_values = svc._compute_dtt(tree, data, ordered_names)
        for v in dtt_values:
            assert v >= -1e-10

    def test_zero_height_tree(self):
        """Star tree with zero branch lengths should return fallback."""
        args = _make_args()
        svc = Dtt(args)
        tree = Phylo.read(StringIO("(A:0,B:0,C:0);"), "newick")
        data = np.array([[1.0], [2.0], [3.0]])
        times, dtt_values = svc._compute_dtt(tree, data, ["A", "B", "C"])
        assert times == [0.0, 1.0]
        assert dtt_values == [1.0, 0.0]

    def test_compute_dtt_fast_path_does_not_call_distance(self, monkeypatch):
        args = _make_args()
        svc = Dtt(args)
        newick = "((A:1,B:1):1,C:2):0.5;"
        ordered_names = ["A", "B", "C"]
        data = np.array([[0.0], [2.0], [4.0]])

        expected = svc._compute_dtt(
            Phylo.read(StringIO(newick), "newick"), data, ordered_names
        )
        tree = Phylo.read(StringIO(newick), "newick")

        def fail_distance(*args, **kwargs):
            raise AssertionError("distance fallback should not be called")

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("get_terminals should not be called")

        monkeypatch.setattr(tree, "distance", fail_distance)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)
        observed = svc._compute_dtt(tree, data, ordered_names)

        np.testing.assert_allclose(observed[0], expected[0])
        np.testing.assert_allclose(observed[1], expected[1])

    def test_prepare_dtt_context_uses_direct_standard_tree_setup(self, monkeypatch):
        args = _make_args()
        svc = Dtt(args)
        newick = "((A:1,B:1):3,(C:2,(D:1,E:1):1):2);"
        ordered_names = ["A", "B", "C", "D", "E"]
        data = np.arange(5.0).reshape(-1, 1)

        expected = svc._compute_dtt(
            Phylo.read(StringIO(newick), "newick"),
            data,
            ordered_names,
        )
        tree = Phylo.read(StringIO(newick), "newick")

        def fail_generic_setup(*_args, **_kwargs):
            raise AssertionError("standard DTT setup should use direct traversal")

        monkeypatch.setattr(type(tree), "depths", fail_generic_setup)
        monkeypatch.setattr(TreeMixin, "find_clades", fail_generic_setup)
        monkeypatch.setattr("Bio.Phylo.Newick.Clade.is_terminal", fail_generic_setup)

        observed = svc._compute_dtt(tree, data, ordered_names)

        np.testing.assert_allclose(observed[0], expected[0])
        np.testing.assert_allclose(observed[1], expected[1])

    def test_compute_dtt_lineage_means_avoid_numpy_dispatch(self, monkeypatch):
        args = _make_args()
        svc = Dtt(args)
        newick = "((A:1,B:1):3,(C:2,(D:1,E:1):1):2);"
        ordered_names = ["A", "B", "C", "D", "E"]
        data = np.arange(5.0).reshape(-1, 1)

        expected = svc._compute_dtt(
            Phylo.read(StringIO(newick), "newick"),
            data,
            ordered_names,
        )

        monkeypatch.setattr(
            dtt_module.np,
            "mean",
            lambda *_args, **_kwargs: pytest.fail(
                "DTT lineage means should use Python list reductions"
            ),
        )
        observed = svc._compute_dtt(
            Phylo.read(StringIO(newick), "newick"),
            data,
            ordered_names,
        )

        np.testing.assert_allclose(observed[0], expected[0])
        np.testing.assert_allclose(observed[1], expected[1])

    def test_prepare_dtt_context_matches_legacy_lineage_scan(self):
        args = _make_args()
        svc = Dtt(args)
        tree = Phylo.read(
            StringIO("((A:1,B:1):3,(C:2,(D:1,E:1):1):2);"),
            "newick",
        )
        ordered_names = ["A", "B", "C", "D", "E"]

        context = svc._prepare_dtt_context(tree, ordered_names)

        depths = tree.depths()
        root = tree.root
        root_depth = depths[root]
        tree_height = max(
            depths[terminal] - root_depth
            for terminal in tree.get_terminals()
        )
        all_clades = list(tree.find_clades(order="preorder"))
        clade_tip_indices = {}
        name_to_idx = {name: idx for idx, name in enumerate(ordered_names)}
        for clade in reversed(all_clades):
            if clade.is_terminal():
                clade_tip_indices[id(clade)] = (name_to_idx[clade.name],)
            else:
                clade_tip_indices[id(clade)] = tuple(
                    idx
                    for child in clade.clades
                    for idx in clade_tip_indices[id(child)]
                )
        node_time = {
            id(clade): (depths[clade] - root_depth) / tree_height
            for clade in all_clades
        }
        parent_map = {}
        for clade in all_clades:
            for child in clade.clades:
                parent_map[id(child)] = clade

        internal_nodes = [clade for clade in all_clades if not clade.is_terminal()]
        internal_nodes.sort(key=lambda clade: node_time[id(clade)])
        expected_lineages = []
        for node in internal_nodes:
            t = node_time[id(node)]
            lineage_ids = []
            for clade in all_clades:
                if clade == root:
                    continue
                parent = parent_map.get(id(clade))
                parent_t = node_time[id(parent)]
                child_t = 1.0 if clade.is_terminal() else node_time[id(clade)]
                if (
                    parent_t <= t + 1e-10
                    and child_t > t + 1e-10
                    and len(clade_tip_indices[id(clade)]) >= 2
                ):
                    lineage_ids.append(id(clade))
            expected_lineages.append(set(lineage_ids))

        assert [
            set(lineage_ids)
            for lineage_ids in context["lineage_clade_ids"]
        ] == expected_lineages


class TestComputeMdi:
    def test_identical_curves_give_zero_mdi(self):
        """MDI should be 0 when observed DTT equals null median."""
        args = _make_args()
        svc = Dtt(args)
        times = [0.0, 0.5, 1.0]
        dtt_values = [1.0, 0.5, 0.0]
        null_median = np.array([1.0, 0.5, 0.0])
        mdi = svc._compute_mdi(times, dtt_values, null_median)
        assert mdi == pytest.approx(0.0, abs=1e-10)

    def test_positive_mdi_for_high_observed(self):
        """MDI should be positive when observed is above null."""
        args = _make_args()
        svc = Dtt(args)
        times = [0.0, 0.5, 1.0]
        dtt_values = [1.0, 0.8, 0.0]
        null_median = np.array([1.0, 0.5, 0.0])
        mdi = svc._compute_mdi(times, dtt_values, null_median)
        assert mdi > 0

    def test_negative_mdi_for_low_observed(self):
        """MDI should be negative when observed is below null."""
        args = _make_args()
        svc = Dtt(args)
        times = [0.0, 0.5, 1.0]
        dtt_values = [1.0, 0.2, 0.0]
        null_median = np.array([1.0, 0.5, 0.0])
        mdi = svc._compute_mdi(times, dtt_values, null_median)
        assert mdi < 0


class TestRun:
    def test_tips_to_prune_for_traits_ordered_all_shared_skips_membership_scan(
        self, monkeypatch
    ):
        tree_tips = ["a", "b", "c"]
        traits = {"a": [1.0], "b": [2.0], "c": [3.0]}

        def fail_iter(*_args, **_kwargs):
            raise AssertionError(
                "ordered all-shared traits should not build a set"
            )

        monkeypatch.setattr("builtins.set", fail_iter)

        assert Dtt._tips_to_prune_for_traits(tree_tips, traits) == []

    def test_tips_to_prune_for_traits_returns_ordered_tail_without_set(
        self, monkeypatch
    ):
        tree_tips = ["a", "b", "c", "d"]
        traits = {"a": [1.0], "b": [2.0], "c": [3.0]}

        def fail_iter(*_args, **_kwargs):
            raise AssertionError("ordered trait prefix should not build a set")

        monkeypatch.setattr("builtins.set", fail_iter)

        assert Dtt._tips_to_prune_for_traits(tree_tips, traits) == ["d"]

    def test_tips_to_prune_for_traits_preserves_interleaved_pruning(self):
        tree_tips = ["a", "d", "b", "c"]
        traits = {"a": [1.0], "b": [2.0], "c": [3.0]}

        assert Dtt._tips_to_prune_for_traits(tree_tips, traits) == ["d"]

    def test_run_all_tips_present_uses_read_only_tree_without_copy_or_prune(
        self, mocker
    ):
        args = _make_args(json=True)
        svc = Dtt(args)
        tree = Phylo.read(StringIO("((a:1,b:1):1,c:2);"), "newick")
        mocked_read = mocker.patch.object(
            Dtt,
            "read_tree_file_unmodified",
            return_value=tree,
        )
        mocker.patch.object(
            Dtt,
            "read_tree_file",
            side_effect=AssertionError("run should not copy cached trees"),
        )
        mocked_validate = mocker.patch.object(Dtt, "validate_tree")
        mocker.patch.object(
            Dtt,
            "get_tip_names_from_tree",
            return_value=["a", "b", "c"],
        )
        mocker.patch(
            "phykit.services.tree.dtt.parse_multi_trait_file",
            return_value=(
                ["x", "y"],
                {
                    "a": [1.0, 2.0],
                    "b": [2.0, 3.0],
                    "c": [3.0, 4.0],
                },
            ),
        )
        mocker.patch.object(
            Dtt,
            "_fast_copy",
            side_effect=AssertionError("all-shared tree should not be copied"),
        )
        mocker.patch.object(
            Dtt,
            "prune_tree_using_taxa_list",
            side_effect=AssertionError("all-shared tree should not be pruned"),
        )
        captured = {}

        def compute_dtt(tree_arg, data_arg, ordered_names_arg, *_args):
            captured["tree"] = tree_arg
            captured["data"] = data_arg
            captured["ordered_names"] = ordered_names_arg
            return [0.0, 1.0], [1.0, 0.0]

        mocker.patch.object(
            Dtt,
            "_compute_dtt",
            side_effect=compute_dtt,
        )
        mocker.patch("phykit.services.tree.dtt.print_json")

        svc.run()

        mocked_read.assert_called_once_with()
        mocked_validate.assert_called_once_with(
            tree,
            min_tips=3,
            require_branch_lengths=True,
            context="disparity through time",
        )
        assert captured["tree"] is tree
        assert captured["ordered_names"] == ["a", "b", "c"]
        assert captured["data"].tolist() == [[1.0, 2.0], [2.0, 3.0], [3.0, 4.0]]

    def test_run_selected_trait_column_builds_one_column_matrix(self, mocker):
        args = _make_args(trait="y", json=True)
        svc = Dtt(args)
        tree = Phylo.read(StringIO("((a:1,b:1):1,c:2);"), "newick")
        mocker.patch.object(Dtt, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(Dtt, "validate_tree")
        mocker.patch.object(
            Dtt,
            "get_tip_names_from_tree",
            return_value=["a", "b", "c"],
        )
        mocker.patch(
            "phykit.services.tree.dtt.parse_multi_trait_file",
            return_value=(
                ["x", "y"],
                {
                    "a": [1.0, 2.0],
                    "b": [2.0, 3.0],
                    "c": [3.0, 4.0],
                },
            ),
        )
        captured = {}

        def compute_dtt(_tree, data_arg, ordered_names_arg, *_args):
            captured["data"] = data_arg
            captured["ordered_names"] = ordered_names_arg
            return [0.0, 1.0], [1.0, 0.0]

        mocker.patch.object(Dtt, "_compute_dtt", side_effect=compute_dtt)
        mocker.patch("phykit.services.tree.dtt.print_json")

        svc.run()

        assert captured["ordered_names"] == ["a", "b", "c"]
        assert captured["data"].tolist() == [[2.0], [3.0], [4.0]]

    def test_run_reuses_initial_tip_names_for_prune_setup(self, mocker):
        args = _make_args(trait="body_mass", json=True)
        svc = Dtt(args)
        tip_name_spy = mocker.spy(svc, "get_tip_names_from_tree")
        mocker.patch.object(
            Dtt,
            "_compute_dtt",
            return_value=([0.0, 1.0], [1.0, 0.0]),
        )
        mocker.patch("phykit.services.tree.dtt.print_json")

        svc.run()

        assert tip_name_spy.call_count == 1

    def test_run_missing_trait_taxa_copies_before_pruning(self, monkeypatch):
        args = _make_args(json=True)
        svc = Dtt(args)
        tree = Phylo.read(StringIO("((a:1,b:1):1,(c:1,d:1):1);"), "newick")
        original_fast_copy = svc._fast_copy
        copied_trees = []
        captured = {}

        def copy_spy(tree_to_copy):
            copied_tree = original_fast_copy(tree_to_copy)
            copied_trees.append(copied_tree)
            return copied_tree

        def compute_dtt(tree_arg, *_args):
            captured["tree"] = tree_arg
            return [0.0, 1.0], [1.0, 0.0]

        monkeypatch.setattr(svc, "read_tree_file_unmodified", lambda: tree)
        monkeypatch.setattr(svc, "_fast_copy", copy_spy)
        monkeypatch.setattr(
            dtt_module,
            "parse_multi_trait_file",
            lambda *_args, **_kwargs: (
                ["x", "y"],
                {
                    "a": [1.0, 2.0],
                    "b": [2.0, 3.0],
                    "c": [3.0, 4.0],
                },
            ),
        )
        monkeypatch.setattr(svc, "_compute_dtt", compute_dtt)
        monkeypatch.setattr(
            dtt_module,
            "print_json",
            lambda *_args, **_kwargs: None,
        )

        svc.run()

        assert len(copied_trees) == 1
        assert captured["tree"] is copied_trees[0]
        assert {tip.name for tip in tree.get_terminals()} == {"a", "b", "c", "d"}
        assert {tip.name for tip in copied_trees[0].get_terminals()} == {
            "a",
            "b",
            "c",
        }

    def test_text_output_single_trait(self, capsys):
        args = _make_args(trait="body_mass")
        svc = Dtt(args)
        svc.run()
        out, _ = capsys.readouterr()
        assert "Disparity Through Time (DTT)" in out
        assert "Index: avg_sq" in out
        assert "N time points:" in out
        assert "Time" in out
        assert "Relative disparity" in out

    def test_print_text_batches_full_report(self, mocker):
        args = _make_args(index="avg_manhattan")
        svc = Dtt(args)
        printed = mocker.patch("builtins.print")

        svc._print_text([0.0, 0.5, 1.0], [1.0, 0.25, 0.0], 0.1234567, 0.04567)

        printed.assert_called_once_with(
            "Disparity Through Time (DTT)\n"
            "Index: avg_manhattan\n"
            "N time points: 3\n"
            "MDI: 0.123457\n"
            "MDI p-value: 0.0457\n"
            "\n"
            "      Time   Relative disparity\n"
            "--------------------------------\n"
            "  0.000000             1.000000\n"
            "  0.500000             0.250000\n"
            "  1.000000             0.000000"
        )

    def test_text_output_all_traits(self, capsys):
        args = _make_args()
        svc = Dtt(args)
        svc.run()
        out, _ = capsys.readouterr()
        assert "Disparity Through Time (DTT)" in out

    def test_json_output(self, capsys):
        args = _make_args(trait="body_mass", json=True)
        svc = Dtt(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "index" in data
        assert "n_time_points" in data
        assert "times" in data
        assert "dtt" in data
        assert data["index"] == "avg_sq"
        assert len(data["times"]) == len(data["dtt"])
        assert data["times"][0] == pytest.approx(0.0)
        assert data["dtt"][0] == pytest.approx(1.0)

    def test_json_output_with_simulations(self, capsys):
        args = _make_args(trait="body_mass", json=True, nsim=10, seed=42)
        svc = Dtt(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "mdi" in data
        assert "mdi_p_value" in data
        assert isinstance(data["mdi"], float)
        assert 0.0 <= data["mdi_p_value"] <= 1.0

    def test_plot_output(self, tmp_path):
        plot_path = str(tmp_path / "dtt_test.png")
        args = _make_args(trait="body_mass", plot_output=plot_path)
        svc = Dtt(args)
        svc.run()
        assert os.path.exists(plot_path)
        assert os.path.getsize(plot_path) > 0

    def test_plot_with_simulations(self, tmp_path):
        plot_path = str(tmp_path / "dtt_sim.png")
        args = _make_args(trait="body_mass", plot_output=plot_path, nsim=10, seed=42)
        svc = Dtt(args)
        svc.run()
        assert os.path.exists(plot_path)
        assert os.path.getsize(plot_path) > 0

    def test_manhattan_index(self, capsys):
        args = _make_args(trait="body_mass", index="avg_manhattan")
        svc = Dtt(args)
        svc.run()
        out, _ = capsys.readouterr()
        assert "Index: avg_manhattan" in out

    def test_invalid_trait_column_raises(self):
        args = _make_args(trait="nonexistent_trait")
        svc = Dtt(args)
        with pytest.raises(SystemExit):
            svc.run()

    def test_too_few_tips_raises(self, tmp_path):
        tree_file = tmp_path / "two_tips.tre"
        tree_file.write_text("(A:1,B:1);\n")
        traits_file = tmp_path / "two_traits.tsv"
        traits_file.write_text("taxon\tval\nA\t1.0\nB\t2.0\n")
        args = _make_args(tree=str(tree_file), traits=str(traits_file))
        svc = Dtt(args)
        with pytest.raises(SystemExit):
            svc.run()


class TestSimulateNull:
    def test_avg_sq_batch_simulation_matches_scalar_loop(self):
        args = _make_args(trait="body_mass", nsim=5, seed=42)
        svc = Dtt(args)
        tree = svc.read_tree_file()

        ordered_names = sorted(
            ["bear", "cat", "dog", "monkey", "raccoon", "sea_lion", "seal", "weasel"]
        )
        data = np.array([
            {"raccoon": 1.04, "bear": 2.39, "sea_lion": 2.30, "seal": 1.88,
             "monkey": 0.60, "cat": 0.56, "weasel": -0.30, "dog": 1.18}[n]
            for n in ordered_names
        ]).reshape(-1, 1)
        obs_times, _ = svc._compute_dtt(tree, data, ordered_names)

        observed = svc._simulate_null(tree, data, ordered_names, obs_times)

        from phykit.services.tree.vcv_utils import build_vcv_matrix
        from phykit.services.tree.dtt import _trapezoid

        rng = np.random.default_rng(args.seed)
        vcv = build_vcv_matrix(tree, ordered_names)
        L = np.linalg.cholesky(vcv)
        context = svc._prepare_dtt_context(tree, ordered_names)
        expected_matrix = np.zeros((args.nsim, len(obs_times)))
        for sim_idx in range(args.nsim):
            z = rng.standard_normal(len(ordered_names))
            sim_data = (L @ z).reshape(-1, 1)
            sim_times, sim_values = svc._compute_dtt(
                tree,
                sim_data,
                ordered_names,
                context=context,
            )
            expected_matrix[sim_idx, :] = np.interp(
                obs_times,
                sim_times,
                sim_values,
            )

        null_median = np.median(expected_matrix, axis=0)
        obs_arr = np.array(
            svc._compute_dtt(tree, data, ordered_names, context=context)[1]
        )
        obs_times_arr = np.array(obs_times)
        expected_mdi = float(_trapezoid(
            np.interp(obs_times_arr, obs_times, obs_arr) - null_median,
            obs_times_arr,
        ))
        expected_sim_mdis = np.zeros(args.nsim)
        for sim_idx in range(args.nsim):
            expected_sim_mdis[sim_idx] = float(_trapezoid(
                expected_matrix[sim_idx, :] - null_median,
                obs_times_arr,
            ))
        expected_p = float(
            np.mean(np.abs(expected_sim_mdis) >= np.abs(expected_mdi))
        )

        np.testing.assert_allclose(observed[0], expected_matrix)
        assert observed[1] == pytest.approx(expected_mdi)
        assert observed[2] == pytest.approx(expected_p)

    def test_avg_sq_batch_simulation_reuses_compute_dtt(self, mocker):
        args = _make_args(trait="body_mass", nsim=5, seed=42)
        svc = Dtt(args)
        tree = svc.read_tree_file()

        ordered_names = sorted(
            ["bear", "cat", "dog", "monkey", "raccoon", "sea_lion", "seal", "weasel"]
        )
        data = np.array([
            {"raccoon": 1.04, "bear": 2.39, "sea_lion": 2.30, "seal": 1.88,
             "monkey": 0.60, "cat": 0.56, "weasel": -0.30, "dog": 1.18}[n]
            for n in ordered_names
        ]).reshape(-1, 1)
        obs_times, _ = svc._compute_dtt(tree, data, ordered_names)
        compute_spy = mocker.spy(svc, "_compute_dtt")

        svc._simulate_null(tree, data, ordered_names, obs_times)

        assert compute_spy.call_count == 1

    def test_avg_sq_batch_simulation_reuses_observed_dtt_values(
        self, mocker
    ):
        args = _make_args(trait="body_mass", nsim=5, seed=42)
        svc = Dtt(args)
        tree = svc.read_tree_file()

        ordered_names = sorted(
            ["bear", "cat", "dog", "monkey", "raccoon", "sea_lion", "seal", "weasel"]
        )
        data = np.array([
            {"raccoon": 1.04, "bear": 2.39, "sea_lion": 2.30, "seal": 1.88,
             "monkey": 0.60, "cat": 0.56, "weasel": -0.30, "dog": 1.18}[n]
            for n in ordered_names
        ]).reshape(-1, 1)
        obs_times, obs_values = svc._compute_dtt(tree, data, ordered_names)
        expected = svc._simulate_null(tree, data, ordered_names, obs_times)

        mocker.patch.object(
            svc,
            "_compute_dtt",
            side_effect=AssertionError("observed DTT values should be reused"),
        )

        observed = svc._simulate_null(
            tree,
            data,
            ordered_names,
            obs_times,
            obs_values,
        )

        np.testing.assert_allclose(observed[0], expected[0])
        assert observed[1] == pytest.approx(expected[1])
        assert observed[2] == pytest.approx(expected[2])

    def test_avg_sq_batch_simulation_preallocates_terminal_point(self, monkeypatch):
        args = _make_args(trait="body_mass", nsim=5, seed=42)
        svc = Dtt(args)
        tree = svc.read_tree_file()

        ordered_names = sorted(
            ["bear", "cat", "dog", "monkey", "raccoon", "sea_lion", "seal", "weasel"]
        )
        data = np.array([
            {"raccoon": 1.04, "bear": 2.39, "sea_lion": 2.30, "seal": 1.88,
             "monkey": 0.60, "cat": 0.56, "weasel": -0.30, "dog": 1.18}[n]
            for n in ordered_names
        ]).reshape(-1, 1)
        context = svc._prepare_dtt_context(tree, ordered_names)
        obs_times = list(context["times"]) + [1.0]
        assert obs_times[-1] == pytest.approx(1.0)

        def fail_append(*_args, **_kwargs):
            raise AssertionError("terminal point should be preallocated")

        monkeypatch.setattr(dtt_module.np, "append", fail_append)
        monkeypatch.setattr(dtt_module.np, "column_stack", fail_append)
        monkeypatch.setattr(dtt_module.np, "where", fail_append)
        monkeypatch.setattr(
            svc,
            "_compute_dtt",
            lambda *_args, **_kwargs: (obs_times, [1.0] * len(obs_times)),
        )

        sim_dtt, _, _ = svc._simulate_null(tree, data, ordered_names, obs_times)

        assert sim_dtt.shape == (args.nsim, len(obs_times))

    def test_avg_sq_batch_simulation_skips_interp_for_strict_matching_grid(
        self, monkeypatch
    ):
        args = _make_args(trait="body_mass", nsim=5, seed=42)
        svc = Dtt(args)
        obs_times = [0.0, 0.5, 1.0]
        context = {
            "zero_height": False,
            "times": obs_times,
            "lineage_clade_ids": [[1, 2], []],
            "clade_index_arrays": {
                1: np.asarray([0, 1], dtype=np.intp),
                2: np.asarray([2, 3], dtype=np.intp),
            },
        }
        data = np.arange(8, dtype=float).reshape(4, 2)
        cholesky_factor = np.eye(4)

        def fail_interp(*_args, **_kwargs):
            raise AssertionError("strict matching time grids should skip interpolation")

        monkeypatch.setattr(dtt_module.np, "interp", fail_interp)

        sim_dtt, _, _ = svc._simulate_null_avg_sq_batch(
            cholesky_factor,
            data,
            obs_times,
            context,
        )

        assert sim_dtt.shape == (args.nsim, len(obs_times))

    def test_avg_sq_clade_disparities_use_postorder_aggregates(self):
        args = _make_args(trait="body_mass", nsim=4, seed=42)
        svc = Dtt(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        ordered_names = ["A", "B", "C", "D"]
        context = svc._prepare_dtt_context(tree, ordered_names)
        clade_ids = list(context["clade_index_arrays"].keys())
        clade_id_to_pos = {
            clade_id: idx for idx, clade_id in enumerate(clade_ids)
        }
        sim_data = np.arange(args.nsim * 4 * 2, dtype=float).reshape(args.nsim, 4, 2)

        expected = np.zeros((args.nsim, len(clade_ids)))
        for clade_id, clade_pos in clade_id_to_pos.items():
            indices = context["clade_index_arrays"][clade_id]
            m = len(indices)
            subset = sim_data[:, indices, :]
            expected[:, clade_pos] = (
                m * np.sum(subset * subset, axis=(1, 2))
                - np.sum(np.sum(subset, axis=1) ** 2, axis=1)
            ) / (m * (m - 1) / 2)

        context_without_index_arrays = dict(context)
        context_without_index_arrays["clade_index_arrays"] = {}
        observed = svc._batch_clade_disparities_avg_sq(
            sim_data,
            context_without_index_arrays,
            clade_id_to_pos,
        )

        np.testing.assert_allclose(observed, expected)

    def test_simulation_reuses_prepared_dtt_context(self, mocker):
        args = _make_args(trait="body_mass", nsim=5, seed=42)
        svc = Dtt(args)
        tree = svc.read_tree_file()

        ordered_names = sorted(
            ["bear", "cat", "dog", "monkey", "raccoon", "sea_lion", "seal", "weasel"]
        )
        data = np.array([
            {"raccoon": 1.04, "bear": 2.39, "sea_lion": 2.30, "seal": 1.88,
             "monkey": 0.60, "cat": 0.56, "weasel": -0.30, "dog": 1.18}[n]
            for n in ordered_names
        ]).reshape(-1, 1)
        times, _ = svc._compute_dtt(tree, data, ordered_names)
        context_spy = mocker.spy(svc, "_prepare_dtt_context")

        svc._simulate_null(tree, data, ordered_names, times)

        assert context_spy.call_count == 1

    def test_simulation_returns_correct_shapes(self):
        args = _make_args(trait="body_mass", nsim=5, seed=42)
        svc = Dtt(args)
        tree = svc.read_tree_file()

        ordered_names = sorted(
            ["bear", "cat", "dog", "monkey", "raccoon", "sea_lion", "seal", "weasel"]
        )
        data = np.array([
            {"raccoon": 1.04, "bear": 2.39, "sea_lion": 2.30, "seal": 1.88,
             "monkey": 0.60, "cat": 0.56, "weasel": -0.30, "dog": 1.18}[n]
            for n in ordered_names
        ]).reshape(-1, 1)

        times, dtt_values = svc._compute_dtt(tree, data, ordered_names)
        sim_dtt, mdi, mdi_p = svc._simulate_null(tree, data, ordered_names, times)

        assert sim_dtt.shape[0] == 5
        assert sim_dtt.shape[1] == len(times)
        assert isinstance(mdi, float)
        assert 0.0 <= mdi_p <= 1.0

    def test_seed_reproducibility(self):
        args1 = _make_args(trait="body_mass", nsim=5, seed=123)
        svc1 = Dtt(args1)
        tree1 = svc1.read_tree_file()
        ordered_names = sorted(
            ["bear", "cat", "dog", "monkey", "raccoon", "sea_lion", "seal", "weasel"]
        )
        data = np.array([
            {"raccoon": 1.04, "bear": 2.39, "sea_lion": 2.30, "seal": 1.88,
             "monkey": 0.60, "cat": 0.56, "weasel": -0.30, "dog": 1.18}[n]
            for n in ordered_names
        ]).reshape(-1, 1)
        times1, _ = svc1._compute_dtt(tree1, data, ordered_names)
        sim1, mdi1, _ = svc1._simulate_null(tree1, data, ordered_names, times1)

        args2 = _make_args(trait="body_mass", nsim=5, seed=123)
        svc2 = Dtt(args2)
        tree2 = svc2.read_tree_file()
        times2, _ = svc2._compute_dtt(tree2, data, ordered_names)
        sim2, mdi2, _ = svc2._simulate_null(tree2, data, ordered_names, times2)

        np.testing.assert_array_almost_equal(sim1, sim2)
        assert mdi1 == pytest.approx(mdi2)
