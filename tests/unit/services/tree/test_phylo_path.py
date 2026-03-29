from argparse import Namespace
import json

import pytest

from phykit.services.tree.phylo_path import PhyloPath


TREE = "tests/sample_files/tree_simple.tre"
TRAITS = "tests/sample_files/tree_simple_multi_traits.tsv"
MODELS = "tests/sample_files/phylo_path_models.txt"


def _make_args(**overrides):
    defaults = dict(
        tree=TREE, traits=TRAITS, models=MODELS,
        best_only=False, plot_output=None, csv=None, json=False,
        fig_width=None, fig_height=None, dpi=150, no_title=False,
        title=None, legend_position=None, ylabel_fontsize=None,
        xlabel_fontsize=None, title_fontsize=None, axis_fontsize=None,
        colors=None, ladderize=False, cladogram=False, circular=False,
        color_file=None,
    )
    defaults.update(overrides)
    return Namespace(**defaults)


class TestPhyloPath:
    def test_basic_run(self, capsys):
        """Runs without error and prints model comparison."""
        svc = PhyloPath(_make_args())
        svc.run()
        captured = capsys.readouterr()
        assert "Phylogenetic Path Analysis" in captured.out
        assert "Model comparison" in captured.out
        assert "Path coefficients" in captured.out

    def test_three_models_ranked(self, capsys):
        """All three models appear in output."""
        svc = PhyloPath(_make_args())
        svc.run()
        captured = capsys.readouterr()
        assert "direct" in captured.out
        assert "mediated" in captured.out
        assert "full" in captured.out

    def test_best_model_identified(self, capsys):
        """Best model is identified."""
        svc = PhyloPath(_make_args())
        svc.run()
        captured = capsys.readouterr()
        assert "Best model:" in captured.out

    def test_model_averaging_default(self, capsys):
        """Default output uses model averaging."""
        svc = PhyloPath(_make_args())
        svc.run()
        captured = capsys.readouterr()
        assert "model-averaged" in captured.out

    def test_best_only_flag(self, capsys):
        """--best-only shows only best model coefficients."""
        svc = PhyloPath(_make_args(best_only=True))
        svc.run()
        captured = capsys.readouterr()
        assert "best model:" in captured.out

    def test_json_structure(self, capsys):
        """JSON output has correct structure."""
        svc = PhyloPath(_make_args(json=True))
        svc.run()
        captured = capsys.readouterr()
        for line in captured.out.strip().split("\n"):
            if line.startswith("{"):
                payload = json.loads(line)
                break
        assert payload["n_taxa"] == 8
        assert payload["n_models"] == 3
        assert "model_comparison" in payload
        assert "path_coefficients" in payload
        assert len(payload["model_comparison"]) == 3

    def test_ciccs_are_finite(self, capsys):
        """CICc values should be finite numbers."""
        svc = PhyloPath(_make_args(json=True))
        svc.run()
        captured = capsys.readouterr()
        for line in captured.out.strip().split("\n"):
            if line.startswith("{"):
                payload = json.loads(line)
                break
        for m in payload["model_comparison"]:
            assert m["CICc"] > 0

    def test_weights_sum_to_one(self, capsys):
        """Model weights should sum to ~1."""
        svc = PhyloPath(_make_args(json=True))
        svc.run()
        captured = capsys.readouterr()
        for line in captured.out.strip().split("\n"):
            if line.startswith("{"):
                payload = json.loads(line)
                break
        total_w = sum(m["weight"] for m in payload["model_comparison"])
        assert total_w == pytest.approx(1.0, abs=0.01)

    def test_path_coefficients_present(self, capsys):
        """Path coefficients are in output."""
        svc = PhyloPath(_make_args(json=True))
        svc.run()
        captured = capsys.readouterr()
        for line in captured.out.strip().split("\n"):
            if line.startswith("{"):
                payload = json.loads(line)
                break
        coefs = payload["path_coefficients"]
        assert "body_mass -> brain_size" in coefs
        # body_mass->brain_size should be positive and strong
        assert coefs["body_mass -> brain_size"]["coef"] > 0.5

    def test_csv_output(self, tmp_path):
        """CSV file is created."""
        csv_path = str(tmp_path / "path.csv")
        svc = PhyloPath(_make_args(csv=csv_path))
        svc.run()
        with open(csv_path) as f:
            content = f.read()
        assert "model,k,q" in content
        assert "path,coef,se" in content

    def test_plot_output(self, tmp_path):
        """Plot file is created."""
        plot_path = str(tmp_path / "dag.png")
        svc = PhyloPath(_make_args(plot_output=plot_path))
        svc.run()
        assert (tmp_path / "dag.png").exists()

    def test_missing_models_file_raises(self):
        """Nonexistent models file raises error."""
        svc = PhyloPath(_make_args(models="nonexistent.txt"))
        with pytest.raises(SystemExit) as exc_info:
            svc.run()
        assert exc_info.value.code == 2

    def test_single_model_raises(self, tmp_path):
        """Only 1 model should raise error (need >= 2)."""
        m = tmp_path / "one_model.txt"
        m.write_text("only: body_mass->brain_size\n")
        svc = PhyloPath(_make_args(models=str(m)))
        with pytest.raises(SystemExit) as exc_info:
            svc.run()
        assert exc_info.value.code == 2
