import json
import os
from argparse import Namespace
from pathlib import Path

import numpy as np
import pytest

from phykit.services.alignment.identity_matrix import IdentityMatrix


def _write_alignment(path, seqs):
    """Write a dict of {name: seq} as FASTA."""
    with open(path, "w") as f:
        for name, seq in seqs.items():
            f.write(f">{name}\n{seq}\n")


def _make_args(alignment, output, **kwargs):
    defaults = dict(
        alignment=str(alignment),
        output=str(output),
        metric="identity",
        tree=None,
        sort="cluster",
        partition=None,
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


class TestIdentityMatrixUnit:
    def test_identity_self_is_one(self, tmp_path):
        """Diagonal entries should all be 1.0."""
        aln_path = tmp_path / "aln.fa"
        out_path = tmp_path / "out.png"
        _write_alignment(aln_path, {
            "taxon_A": "ACGTACGT",
            "taxon_B": "ACGTTTTT",
            "taxon_C": "TTTTACGT",
        })
        args = _make_args(aln_path, out_path)
        service = IdentityMatrix(args)
        sequences, taxa_names, aln_length = service._parse_alignment()
        matrix = service._compute_identity_matrix(sequences, taxa_names)
        for i in range(len(taxa_names)):
            assert matrix[i, i] == 1.0

    def test_identity_symmetric(self, tmp_path):
        """Matrix[i,j] should equal matrix[j,i]."""
        aln_path = tmp_path / "aln.fa"
        out_path = tmp_path / "out.png"
        _write_alignment(aln_path, {
            "taxon_A": "ACGTACGT",
            "taxon_B": "ACGTTTTT",
            "taxon_C": "TTTTACGT",
        })
        args = _make_args(aln_path, out_path)
        service = IdentityMatrix(args)
        sequences, taxa_names, aln_length = service._parse_alignment()
        matrix = service._compute_identity_matrix(sequences, taxa_names)
        n = len(taxa_names)
        for i in range(n):
            for j in range(n):
                assert matrix[i, j] == matrix[j, i]

    def test_identical_sequences(self, tmp_path):
        """Two identical sequences should have identity 1.0."""
        aln_path = tmp_path / "aln.fa"
        out_path = tmp_path / "out.png"
        _write_alignment(aln_path, {
            "taxon_A": "ACGTACGT",
            "taxon_B": "ACGTACGT",
        })
        args = _make_args(aln_path, out_path)
        service = IdentityMatrix(args)
        sequences, taxa_names, aln_length = service._parse_alignment()
        matrix = service._compute_identity_matrix(sequences, taxa_names)
        assert matrix[0, 1] == 1.0

    def test_completely_different(self, tmp_path):
        """No matching non-ambiguous positions should give identity 0.0."""
        aln_path = tmp_path / "aln.fa"
        out_path = tmp_path / "out.png"
        _write_alignment(aln_path, {
            "taxon_A": "AAAA",
            "taxon_B": "TTTT",
        })
        args = _make_args(aln_path, out_path)
        service = IdentityMatrix(args)
        sequences, taxa_names, aln_length = service._parse_alignment()
        matrix = service._compute_identity_matrix(sequences, taxa_names)
        assert matrix[0, 1] == 0.0

    def test_p_distance_metric(self, tmp_path):
        """p-distance should be 1 - identity."""
        aln_path = tmp_path / "aln.fa"
        out_path = tmp_path / "out.png"
        _write_alignment(aln_path, {
            "taxon_A": "ACGT",
            "taxon_B": "ACTT",
        })
        args = _make_args(aln_path, out_path, metric="p-distance")
        service = IdentityMatrix(args)
        sequences, taxa_names, aln_length = service._parse_alignment()
        identity_matrix = service._compute_identity_matrix(sequences, taxa_names)
        # identity should be 3/4 = 0.75
        assert abs(identity_matrix[0, 1] - 0.75) < 1e-9
        # p-distance = 1 - 0.75 = 0.25
        p_dist_matrix = 1.0 - identity_matrix
        assert abs(p_dist_matrix[0, 1] - 0.25) < 1e-9

    def test_gaps_excluded(self, tmp_path):
        """Gaps and ambiguous chars should be skipped in comparison."""
        aln_path = tmp_path / "aln.fa"
        out_path = tmp_path / "out.png"
        # Positions: A=A (match), C=C (match), -=G (gap in A, skip), T=T (match)
        # Valid positions: 3, matches: 3, identity = 1.0
        _write_alignment(aln_path, {
            "taxon_A": "AC-T",
            "taxon_B": "ACGT",
        })
        args = _make_args(aln_path, out_path)
        service = IdentityMatrix(args)
        sequences, taxa_names, aln_length = service._parse_alignment()
        matrix = service._compute_identity_matrix(sequences, taxa_names)
        # Position 3 has a gap in taxon_A, so only 3 valid positions, all match
        assert matrix[0, 1] == 1.0

    def test_sort_alpha(self, tmp_path):
        """Alphabetical ordering should sort taxa by name."""
        aln_path = tmp_path / "aln.fa"
        out_path = tmp_path / "out.png"
        _write_alignment(aln_path, {
            "Zebra": "ACGT",
            "Apple": "ACGT",
            "Mango": "ACGT",
        })
        args = _make_args(aln_path, out_path, sort="alpha")
        service = IdentityMatrix(args)
        sequences, taxa_names, aln_length = service._parse_alignment()
        matrix = service._compute_identity_matrix(sequences, taxa_names)
        order, ordered_labels = service._determine_order(matrix, taxa_names, len(taxa_names))
        assert ordered_labels == ["Apple", "Mango", "Zebra"]

    def test_sort_tree(self, tmp_path):
        """Tree-guided ordering should follow tree tip order."""
        aln_path = tmp_path / "aln.fa"
        out_path = tmp_path / "out.png"
        tree_path = tmp_path / "tree.nwk"
        _write_alignment(aln_path, {
            "A": "ACGT",
            "B": "ACTT",
            "C": "TTTT",
        })
        with open(tree_path, "w") as f:
            f.write("((C,A),B);")
        args = _make_args(aln_path, out_path, sort="tree", tree=str(tree_path))
        service = IdentityMatrix(args)
        sequences, taxa_names, aln_length = service._parse_alignment()
        matrix = service._compute_identity_matrix(sequences, taxa_names)
        order, ordered_labels = service._determine_order(matrix, taxa_names, len(taxa_names))
        assert ordered_labels == ["C", "A", "B"]

    def test_sort_tree_requires_tree(self, tmp_path):
        """--sort tree without --tree should raise SystemExit."""
        aln_path = tmp_path / "aln.fa"
        out_path = tmp_path / "out.png"
        _write_alignment(aln_path, {
            "A": "ACGT",
            "B": "ACTT",
        })
        args = _make_args(aln_path, out_path, sort="tree", tree=None)
        service = IdentityMatrix(args)
        sequences, taxa_names, aln_length = service._parse_alignment()
        matrix = service._compute_identity_matrix(sequences, taxa_names)
        with pytest.raises(SystemExit):
            service._determine_order(matrix, taxa_names, len(taxa_names))

    def test_creates_png(self, tmp_path):
        """The run method should create a non-empty output file."""
        aln_path = tmp_path / "aln.fa"
        out_path = tmp_path / "out.png"
        _write_alignment(aln_path, {
            "taxon_A": "ACGTACGT",
            "taxon_B": "ACGTTTTT",
            "taxon_C": "TTTTACGT",
        })
        args = _make_args(aln_path, out_path)
        service = IdentityMatrix(args)
        service.run()
        assert out_path.exists()
        assert out_path.stat().st_size > 0

    def test_json_output(self, tmp_path, capsys):
        """JSON output should contain correct structure and values."""
        aln_path = tmp_path / "aln.fa"
        out_path = tmp_path / "out.png"
        _write_alignment(aln_path, {
            "taxon_A": "ACGT",
            "taxon_B": "ACGT",
            "taxon_C": "TTTT",
        })
        args = _make_args(aln_path, out_path, json=True)
        service = IdentityMatrix(args)
        service.run()
        captured = capsys.readouterr()
        payload = json.loads(captured.out)
        assert payload["n_taxa"] == 3
        assert payload["alignment_length"] == 4
        assert payload["metric"] == "identity"
        assert "mean_identity" in payload
        assert "min_identity" in payload
        assert "max_identity" in payload
        assert "taxa_order" in payload
        assert "output_file" in payload
        assert len(payload["taxa_order"]) == 3
        # taxon_A and taxon_B are identical, so max should be 1.0
        assert payload["max_identity"]["value"] == 1.0
