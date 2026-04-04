"""
Unit tests for disparity through time (DTT).

Computes DTT curves and MDI following Harmon et al. (2003).
Cross-validated against R's geiger::dtt().

See tests/r_validation/validate_dtt.R for the R validation script.
"""
import pytest
import json
import os
import numpy as np
from argparse import Namespace
from io import StringIO
from pathlib import Path

from Bio import Phylo

from phykit.services.tree.dtt import Dtt


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
