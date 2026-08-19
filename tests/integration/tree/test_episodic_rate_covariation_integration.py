"""End-to-end tests for episodic evolutionary-rate covariation."""

import json
from pathlib import Path
import sys
from unittest.mock import patch

from Bio import Phylo
import pytest

from phykit.phykit import Phykit


SAMPLE_FILES = Path(__file__).parents[2] / "sample_files"


def _run_cli(arguments):
    with patch.object(sys, "argv", ["phykit", *arguments]):
        Phykit()


@pytest.mark.integration
def test_episodic_rate_covariation_reports_corrected_clade_scans(capsys):
    _run_cli(
        [
            "episodic_rate_covariation",
            str(SAMPLE_FILES / "tree_simple.tre"),
            str(SAMPLE_FILES / "tree_simple_other_topology.tre"),
            "--reference",
            str(SAMPLE_FILES / "tree_simple_2.tre"),
            "--permutations",
            "19",
            "--min-edges",
            "2",
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert payload["method"] == "episodic_rate_covariation"
    assert payload["experimental"] is True
    assert payload["permutations"] == 19
    assert payload["candidate_count"] > 0
    assert payload["scan_edge_count"] <= payload["retained_edge_count"]
    assert len(payload["clade_scans"]) == payload["candidate_count"]
    assert all(
        row["fwer_p_value"] >= row["unadjusted_p_value"]
        for row in payload["clade_scans"]
    )
    assert any(
        row["local_contribution"] is not None for row in payload["branches"]
    )


@pytest.mark.integration
def test_erc_scan_alias_writes_table_and_annotated_tree(tmp_path, capsys):
    table_path = tmp_path / "episodic_clades.tsv"
    tree_path = tmp_path / "episodic_reference.tre"
    _run_cli(
        [
            "erc_scan",
            str(SAMPLE_FILES / "tree_simple.tre"),
            str(SAMPLE_FILES / "tree_simple_other_topology.tre"),
            "-r",
            str(SAMPLE_FILES / "tree_simple_2.tre"),
            "--permutations",
            "9",
            "--min-edges",
            "2",
            "--alpha",
            "1",
            "--output",
            str(table_path),
            "--annotated-tree",
            str(tree_path),
        ]
    )

    output = capsys.readouterr().out
    assert "Episodic Rate Covariation (experimental)" in output
    assert table_path.read_text().startswith("clade\ttaxon_count\troot_depth")
    annotated = Phylo.read(tree_path, "newick")
    assert any(clade.comment for clade in annotated.get_nonterminals())
