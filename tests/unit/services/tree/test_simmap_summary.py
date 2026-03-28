from argparse import Namespace
import json

import pytest
import numpy as np

from phykit.services.tree.simmap_summary import SimmapSummary


TREE = "tests/sample_files/tree_simple.tre"
TRAITS = "tests/sample_files/tree_simple_discrete_traits.tsv"


def _make_args(**overrides):
    defaults = dict(
        tree=TREE, trait_data=TRAITS, trait="diet",
        model="ER", nsim=20, seed=42, plot=None, csv=None, json=False,
        fig_width=None, fig_height=None, dpi=150, no_title=False,
        title=None, legend_position=None, ylabel_fontsize=None,
        xlabel_fontsize=None, title_fontsize=None, axis_fontsize=None,
        colors=None, ladderize=False, cladogram=False, circular=False,
        color_file=None,
    )
    defaults.update(overrides)
    return Namespace(**defaults)


class TestSimmapSummary:
    def test_basic_run(self, capsys):
        """Runs without error and prints summary."""
        svc = SimmapSummary(_make_args())
        svc.run()
        captured = capsys.readouterr()
        assert "SIMMAP Summary" in captured.out
        assert "Q matrix:" in captured.out
        assert "Per-branch dwelling" in captured.out
        assert "Node posterior" in captured.out

    def test_log_likelihood_matches_r(self, capsys):
        """Log-likelihood matches phytools::fitMk (ER) closely."""
        svc = SimmapSummary(_make_args())
        svc.run()
        captured = capsys.readouterr()
        # R gives -8.7889; PhyKIT should be close
        assert "-8.78" in captured.out or "-8.79" in captured.out

    def test_all_states_in_output(self, capsys):
        """All states appear in output."""
        svc = SimmapSummary(_make_args())
        svc.run()
        captured = capsys.readouterr()
        assert "carnivore" in captured.out
        assert "herbivore" in captured.out
        assert "omnivore" in captured.out

    def test_per_branch_proportions_sum_to_one(self, capsys):
        """Dwelling proportions on each branch should sum to ~1."""
        svc = SimmapSummary(_make_args(json=True))
        svc.run()
        captured = capsys.readouterr()
        for line in captured.out.strip().split("\n"):
            if line.startswith("{"):
                payload = json.loads(line)
                break
        for branch in payload["branches"]:
            total = sum(branch["dwelling_proportions"].values())
            assert total == pytest.approx(1.0, abs=0.01), (
                f"Branch {branch['branch']} proportions sum to {total}"
            )

    def test_node_posteriors_sum_to_one(self, capsys):
        """Posteriors at each node should sum to 1."""
        svc = SimmapSummary(_make_args(json=True))
        svc.run()
        captured = capsys.readouterr()
        for line in captured.out.strip().split("\n"):
            if line.startswith("{"):
                payload = json.loads(line)
                break
        for node in payload["node_posteriors"]:
            total = sum(node["posteriors"].values())
            assert total == pytest.approx(1.0, abs=0.01)

    def test_json_output_structure(self, capsys):
        """JSON has correct top-level keys."""
        svc = SimmapSummary(_make_args(json=True))
        svc.run()
        captured = capsys.readouterr()
        for line in captured.out.strip().split("\n"):
            if line.startswith("{"):
                payload = json.loads(line)
                break
        assert payload["model"] == "ER"
        assert payload["nsim"] == 20
        assert "q_matrix" in payload
        assert "branches" in payload
        assert "node_posteriors" in payload
        assert "total_expected_transitions" in payload
        assert len(payload["states"]) == 3

    def test_csv_output(self, tmp_path):
        """CSV file is created with expected content."""
        csv_path = str(tmp_path / "simmap.csv")
        svc = SimmapSummary(_make_args(csv=csv_path))
        svc.run()
        with open(csv_path) as f:
            content = f.read()
        assert "branch,branch_length" in content
        assert "prop_carnivore" in content
        assert "posterior_carnivore" in content

    def test_plot_output(self, tmp_path):
        """Plot file is created."""
        plot_path = str(tmp_path / "simmap_pie.png")
        svc = SimmapSummary(_make_args(plot=plot_path))
        svc.run()
        assert (tmp_path / "simmap_pie.png").exists()

    def test_reproducible_with_seed(self, capsys):
        """Same seed gives same results."""
        svc1 = SimmapSummary(_make_args(json=True, seed=123))
        svc1.run()
        out1 = capsys.readouterr().out

        svc2 = SimmapSummary(_make_args(json=True, seed=123))
        svc2.run()
        out2 = capsys.readouterr().out

        j1 = j2 = None
        for line in out1.strip().split("\n"):
            if line.startswith("{"):
                j1 = json.loads(line)
        for line in out2.strip().split("\n"):
            if line.startswith("{"):
                j2 = json.loads(line)
        assert j1["total_expected_transitions"] == j2["total_expected_transitions"]
        assert j1["branches"][0] == j2["branches"][0]

    def test_ard_model(self, capsys):
        """ARD model runs without error."""
        svc = SimmapSummary(_make_args(model="ARD", nsim=10))
        svc.run()
        captured = capsys.readouterr()
        assert "SIMMAP Summary" in captured.out

    def test_expected_transitions_positive(self, capsys):
        """Total expected transitions should be positive."""
        svc = SimmapSummary(_make_args(json=True))
        svc.run()
        captured = capsys.readouterr()
        for line in captured.out.strip().split("\n"):
            if line.startswith("{"):
                payload = json.loads(line)
                break
        assert payload["total_expected_transitions"] > 0

    def test_tip_branches_have_correct_dominant_state(self, capsys):
        """Terminal branches should be dominated by the observed tip state."""
        svc = SimmapSummary(_make_args(json=True, nsim=50))
        svc.run()
        captured = capsys.readouterr()
        for line in captured.out.strip().split("\n"):
            if line.startswith("{"):
                payload = json.loads(line)
                break
        # raccoon is carnivore - its branch should be mostly carnivore
        for b in payload["branches"]:
            if b["branch"] == "raccoon":
                assert b["dwelling_proportions"]["carnivore"] > 0.4
                break
