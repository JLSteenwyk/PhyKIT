from argparse import Namespace
import json

import pytest
import numpy as np

from phykit.services.tree.phylo_anova import PhyloAnova


TREE = "tests/sample_files/tree_simple.tre"
TRAITS_SINGLE = "tests/sample_files/tree_simple_anova_single.tsv"
TRAITS_MULTI = "tests/sample_files/tree_simple_anova_traits.tsv"


def _make_args(**overrides):
    defaults = dict(
        tree=TREE,
        traits=TRAITS_SINGLE,
        group_column=None,
        method="auto",
        permutations=99,
        pairwise=False,
        plot_output=None,
        plot_type="auto",
        seed=42,
        json=False,
        fig_width=None,
        fig_height=None,
        dpi=150,
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
    defaults.update(overrides)
    return Namespace(**defaults)


class TestPhyloAnova:
    def test_anova_deterministic_values(self, capsys):
        """Deterministic SS/MS/F match R validation values."""
        args = _make_args(permutations=99, seed=42)
        svc = PhyloAnova(args)
        svc.run()
        captured = capsys.readouterr()
        # SS_model = 0.0382, SS_resid = 0.2691, F = 0.3548
        assert "0.0382" in captured.out
        assert "0.2691" in captured.out
        assert "0.3548" in captured.out

    def test_anova_json_output(self, capsys):
        """JSON output has correct structure and deterministic values."""
        args = _make_args(json=True, permutations=99, seed=42)
        svc = PhyloAnova(args)
        svc.run()
        captured = capsys.readouterr()
        # Find the JSON line
        for line in captured.out.strip().split("\n"):
            if line.startswith("{"):
                payload = json.loads(line)
                break
        assert payload["method"] == "anova"
        assert payload["n_taxa"] == 8
        assert payload["n_groups"] == 3
        tab = payload["anova_table"]
        assert tab["ss_model"] == pytest.approx(0.038184, abs=1e-4)
        assert tab["ss_resid"] == pytest.approx(0.269069, abs=1e-4)
        assert tab["f_stat"] == pytest.approx(0.354776, abs=1e-3)
        assert tab["df_model"] == 2
        assert tab["df_resid"] == 5

    def test_manova_auto_detection(self, capsys):
        """Multi-trait file auto-detects MANOVA."""
        args = _make_args(traits=TRAITS_MULTI, permutations=99, seed=42)
        svc = PhyloAnova(args)
        svc.run()
        captured = capsys.readouterr()
        assert "MANOVA" in captured.out
        assert "Pillai" in captured.out

    def test_manova_deterministic_values(self, capsys):
        """MANOVA Pillai's trace matches R validation."""
        args = _make_args(
            traits=TRAITS_MULTI, permutations=99, seed=42, json=True,
        )
        svc = PhyloAnova(args)
        svc.run()
        captured = capsys.readouterr()
        for line in captured.out.strip().split("\n"):
            if line.startswith("{"):
                payload = json.loads(line)
                break
        tab = payload["anova_table"]
        assert tab["ss_model"] == pytest.approx(0.049484, abs=1e-4)
        assert tab["ss_resid"] == pytest.approx(0.415777, abs=1e-4)
        assert tab["pillai_trace"] == pytest.approx(0.794354, abs=1e-3)

    def test_pairwise_output(self, capsys):
        """--pairwise produces pairwise comparisons."""
        args = _make_args(pairwise=True, permutations=99, seed=42)
        svc = PhyloAnova(args)
        svc.run()
        captured = capsys.readouterr()
        assert "Pairwise comparisons" in captured.out
        assert "carnivore vs herbivore" in captured.out
        assert "carnivore vs omnivore" in captured.out
        assert "herbivore vs omnivore" in captured.out

    def test_pairwise_json(self, capsys):
        """Pairwise results in JSON output."""
        args = _make_args(pairwise=True, permutations=99, seed=42, json=True)
        svc = PhyloAnova(args)
        svc.run()
        captured = capsys.readouterr()
        for line in captured.out.strip().split("\n"):
            if line.startswith("{"):
                payload = json.loads(line)
                break
        assert "pairwise" in payload
        assert len(payload["pairwise"]) == 3
        groups = {(pw["group1"], pw["group2"]) for pw in payload["pairwise"]}
        assert ("carnivore", "herbivore") in groups

    def test_method_anova_with_multi_trait_raises(self):
        """--method anova with multiple traits raises error."""
        args = _make_args(traits=TRAITS_MULTI, method="anova")
        svc = PhyloAnova(args)
        with pytest.raises(SystemExit) as exc_info:
            svc.run()
        assert exc_info.value.code == 2

    def test_missing_trait_file_raises(self):
        """Nonexistent trait file raises error."""
        args = _make_args(traits="nonexistent.tsv")
        svc = PhyloAnova(args)
        with pytest.raises(SystemExit) as exc_info:
            svc.run()
        assert exc_info.value.code == 2

    def test_bad_group_column_raises(self):
        """Nonexistent group column raises error."""
        args = _make_args(group_column="nonexistent")
        svc = PhyloAnova(args)
        with pytest.raises(SystemExit) as exc_info:
            svc.run()
        assert exc_info.value.code == 2

    def test_plot_boxplot_creates_file(self, tmp_path):
        """Boxplot is created."""
        out = str(tmp_path / "anova_boxplot.png")
        args = _make_args(
            plot_output=out, plot_type="boxplot",
            permutations=99, seed=42,
        )
        svc = PhyloAnova(args)
        svc.run()
        assert (tmp_path / "anova_boxplot.png").exists()

    def test_plot_phylomorphospace_creates_file(self, tmp_path):
        """Phylomorphospace is created for MANOVA."""
        out = str(tmp_path / "manova_morpho.png")
        args = _make_args(
            traits=TRAITS_MULTI,
            plot_output=out, plot_type="phylomorphospace",
            permutations=99, seed=42,
        )
        svc = PhyloAnova(args)
        svc.run()
        assert (tmp_path / "manova_morpho.png").exists()

    def test_reproducible_with_seed(self, capsys):
        """Same seed gives same results."""
        args1 = _make_args(permutations=99, seed=123, json=True)
        svc1 = PhyloAnova(args1)
        svc1.run()
        out1 = capsys.readouterr().out

        args2 = _make_args(permutations=99, seed=123, json=True)
        svc2 = PhyloAnova(args2)
        svc2.run()
        out2 = capsys.readouterr().out

        # Extract JSON lines
        j1 = j2 = None
        for line in out1.strip().split("\n"):
            if line.startswith("{"):
                j1 = json.loads(line)
        for line in out2.strip().split("\n"):
            if line.startswith("{"):
                j2 = json.loads(line)
        assert j1["anova_table"]["p_value"] == j2["anova_table"]["p_value"]
        assert j1["anova_table"]["z_score"] == j2["anova_table"]["z_score"]
