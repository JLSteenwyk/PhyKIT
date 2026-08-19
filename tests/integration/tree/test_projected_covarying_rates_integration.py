"""End-to-end tests for projected evolutionary-rate covariation."""

import json
from pathlib import Path
import sys
from unittest.mock import patch

import pytest

from phykit.phykit import Phykit


SAMPLE_FILES = Path(__file__).parents[2] / "sample_files"


def _run_cli(arguments):
    with patch.object(sys, "argv", ["phykit", *arguments]):
        Phykit()


@pytest.mark.integration
def test_projected_covarying_rates_accepts_discordant_topology(capsys):
    _run_cli(
        [
            "projected_covarying_rates",
            str(SAMPLE_FILES / "tree_simple.tre"),
            str(SAMPLE_FILES / "tree_simple_other_topology.tre"),
            "--reference",
            str(SAMPLE_FILES / "tree_simple_2.tre"),
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert payload["method"] == "projected_covarying_rates"
    assert payload["experimental"] is True
    assert payload["shared_taxon_count"] == 8
    assert payload["reference_edge_count"] == 13
    assert payload["tree_zero_projection"]["nrmse"] == pytest.approx(0.0, abs=1e-12)
    assert payload["tree_one_projection"]["nrmse"] > 0.07
    assert -1.0 <= payload["correlation"] <= 1.0


@pytest.mark.integration
def test_pcover_alias_and_inverse_reference_weighting(capsys):
    _run_cli(
        [
            "pcover",
            str(SAMPLE_FILES / "tree_simple.tre"),
            str(SAMPLE_FILES / "tree_simple_other_topology.tre"),
            "-r",
            str(SAMPLE_FILES / "tree_simple_2.tre"),
            "--weighting",
            "inverse_reference",
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert payload["weighting"] == "inverse_reference"
    assert payload["distance_pairs_used"] == payload["distance_pair_count"]


@pytest.mark.integration
def test_projected_covarying_rates_prunes_to_three_way_taxon_overlap(
    tmp_path,
    capsys,
):
    gene_zero = tmp_path / "gene_zero.tre"
    gene_one = tmp_path / "gene_one.tre"
    reference = tmp_path / "reference.tre"
    gene_zero.write_text(
        "((A:1,B:2):2,((C:1,D:1.5):1,E:2):2);\n",
        encoding="utf-8",
    )
    gene_one.write_text(
        "((A:1,C:1.5):1,(B:2,(D:1,F:2):1):2);\n",
        encoding="utf-8",
    )
    reference.write_text(
        "((A:1,B:1):2,((C:1,D:1):1,(E:1,F:1):1):2);\n",
        encoding="utf-8",
    )

    _run_cli(
        [
            "pcover",
            str(gene_zero),
            str(gene_one),
            "-r",
            str(reference),
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert payload["shared_taxa"] == ["A", "B", "C", "D"]
    assert payload["shared_taxon_count"] == 4


@pytest.mark.integration
def test_projected_covarying_rates_verbose_output_includes_diagnostics(capsys):
    _run_cli(
        [
            "pcover",
            str(SAMPLE_FILES / "tree_simple.tre"),
            str(SAMPLE_FILES / "tree_simple_other_topology.tre"),
            "-r",
            str(SAMPLE_FILES / "tree_simple_2.tre"),
            "--verbose",
        ]
    )

    output = capsys.readouterr().out
    assert "Projected Covarying Evolutionary Rates (experimental)" in output
    assert "Tree one projection NRMSE:" in output
    assert "branch\treference_length\tprojected_zero" in output


@pytest.mark.integration
def test_projected_covarying_rates_writes_plot(tmp_path, capsys):
    plot_path = tmp_path / "projected_rates.png"
    _run_cli(
        [
            "pcover",
            str(SAMPLE_FILES / "tree_simple.tre"),
            str(SAMPLE_FILES / "tree_simple_other_topology.tre"),
            "-r",
            str(SAMPLE_FILES / "tree_simple_2.tre"),
            "--plot",
            "--plot-output",
            str(plot_path),
        ]
    )

    output = capsys.readouterr().out
    assert plot_path.is_file()
    assert plot_path.stat().st_size > 0
    assert f"Saved projected covarying rates plot: {plot_path}" in output
