"""Integration tests for the phylo_gwas command."""

import json
import os
import sys
import tempfile
from pathlib import Path

import pytest
from mock import patch

from phykit.phykit import Phykit


def _write_alignment(path, seqs):
    with open(path, "w") as f:
        for name, seq in seqs.items():
            f.write(f">{name}\n{seq}\n")


def _write_phenotype(path, pheno):
    with open(path, "w") as f:
        for taxon, value in pheno.items():
            f.write(f"{taxon}\t{value}\n")


def _write_tree(path, newick):
    with open(path, "w") as f:
        f.write(newick + "\n")


def _write_partition(path, partitions):
    with open(path, "w") as f:
        for dtype, name, start, end in partitions:
            f.write(f"{dtype}, {name} = {start}-{end}\n")


SEQS = {
    "sp1": "AAGAA",
    "sp2": "AAGCA",
    "sp3": "AAGAA",
    "sp4": "AACAA",
    "sp5": "AGCGA",
    "sp6": "AGCAA",
    "sp7": "AGCGA",
    "sp8": "AGCAA",
}

PHENO = {
    "sp1": "highland",
    "sp2": "highland",
    "sp3": "highland",
    "sp4": "highland",
    "sp5": "lowland",
    "sp6": "lowland",
    "sp7": "lowland",
    "sp8": "lowland",
}

TREE = "((sp1,sp2,sp3,sp4),(sp5,sp6,sp7,sp8));"


@pytest.mark.integration
class TestPhyloGwasIntegration:
    def test_basic_invocation(self, tmp_path):
        aln_path = str(tmp_path / "test.fa")
        pheno_path = str(tmp_path / "pheno.tsv")
        output_path = str(tmp_path / "manhattan.png")

        _write_alignment(aln_path, SEQS)
        _write_phenotype(pheno_path, PHENO)

        testargs = [
            "phykit",
            "phylo_gwas",
            "-a", aln_path,
            "-d", pheno_path,
            "-o", output_path,
        ]
        with patch("builtins.print") as mocked_print:
            with patch.object(sys, "argv", testargs):
                Phykit()

        assert os.path.exists(output_path)
        # Check that text output includes summary
        calls = [str(c) for c in mocked_print.mock_calls]
        call_text = " ".join(calls)
        assert "Phylogenetic GWAS" in call_text

    def test_alias_pgwas(self, tmp_path):
        aln_path = str(tmp_path / "test.fa")
        pheno_path = str(tmp_path / "pheno.tsv")
        output_path = str(tmp_path / "manhattan.png")

        _write_alignment(aln_path, SEQS)
        _write_phenotype(pheno_path, PHENO)

        testargs = [
            "phykit",
            "pgwas",
            "-a", aln_path,
            "-d", pheno_path,
            "-o", output_path,
        ]
        with patch("builtins.print") as mocked_print:
            with patch.object(sys, "argv", testargs):
                Phykit()

        assert os.path.exists(output_path)

    def test_with_tree(self, tmp_path):
        aln_path = str(tmp_path / "test.fa")
        pheno_path = str(tmp_path / "pheno.tsv")
        tree_path = str(tmp_path / "tree.nwk")
        output_path = str(tmp_path / "manhattan.png")

        _write_alignment(aln_path, SEQS)
        _write_phenotype(pheno_path, PHENO)
        _write_tree(tree_path, TREE)

        testargs = [
            "phykit",
            "phylo_gwas",
            "-a", aln_path,
            "-d", pheno_path,
            "-o", output_path,
            "-t", tree_path,
        ]
        with patch("builtins.print") as mocked_print:
            with patch.object(sys, "argv", testargs):
                Phykit()

        assert os.path.exists(output_path)
        calls = [str(c) for c in mocked_print.mock_calls]
        call_text = " ".join(calls)
        # With tree, should report monophyletic/polyphyletic counts
        assert "olyphyletic" in call_text or "onophyletic" in call_text

    def test_with_partition(self, tmp_path):
        aln_path = str(tmp_path / "test.fa")
        pheno_path = str(tmp_path / "pheno.tsv")
        part_path = str(tmp_path / "test.partition")
        output_path = str(tmp_path / "manhattan.png")

        _write_alignment(aln_path, SEQS)
        _write_phenotype(pheno_path, PHENO)
        _write_partition(part_path, [
            ("DNA", "gene1", 1, 3),
            ("DNA", "gene2", 4, 5),
        ])

        testargs = [
            "phykit",
            "phylo_gwas",
            "-a", aln_path,
            "-d", pheno_path,
            "-o", output_path,
            "-p", part_path,
        ]
        with patch("builtins.print"):
            with patch.object(sys, "argv", testargs):
                Phykit()

        assert os.path.exists(output_path)

    def test_json_output(self, tmp_path):
        aln_path = str(tmp_path / "test.fa")
        pheno_path = str(tmp_path / "pheno.tsv")
        output_path = str(tmp_path / "manhattan.png")

        _write_alignment(aln_path, SEQS)
        _write_phenotype(pheno_path, PHENO)

        testargs = [
            "phykit",
            "phylo_gwas",
            "-a", aln_path,
            "-d", pheno_path,
            "-o", output_path,
            "--json",
        ]
        with patch("builtins.print") as mocked_print:
            with patch.object(sys, "argv", testargs):
                Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert "n_taxa" in payload
        assert "alignment_length" in payload
        assert "phenotype_type" in payload
        assert "results" in payload
        assert payload["n_taxa"] == 8
        assert payload["phenotype_type"] == "categorical"

    def test_csv_output(self, tmp_path):
        aln_path = str(tmp_path / "test.fa")
        pheno_path = str(tmp_path / "pheno.tsv")
        output_path = str(tmp_path / "manhattan.png")
        csv_path = str(tmp_path / "results.csv")

        _write_alignment(aln_path, SEQS)
        _write_phenotype(pheno_path, PHENO)

        testargs = [
            "phykit",
            "phylo_gwas",
            "-a", aln_path,
            "-d", pheno_path,
            "-o", output_path,
            "--csv", csv_path,
        ]
        with patch("builtins.print"):
            with patch.object(sys, "argv", testargs):
                Phykit()

        assert os.path.exists(csv_path)
        import csv as csv_mod
        with open(csv_path) as f:
            reader = csv_mod.DictReader(f)
            rows = list(reader)
            assert len(rows) > 0
            assert "position" in reader.fieldnames
            assert "p_value" in reader.fieldnames
