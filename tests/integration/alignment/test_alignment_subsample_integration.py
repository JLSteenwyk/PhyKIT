import json
import os
import pytest
import sys
import tempfile
from mock import patch, call
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)
sample_dir = here.parent.parent.parent / "sample_files"


def _read_lines(path):
    with open(path) as fh:
        return [l.strip() for l in fh if l.strip()]


def _read_fasta(path):
    seqs = {}
    current = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                current = line[1:]
                seqs[current] = ""
            elif current is not None:
                seqs[current] += line
    return seqs


@pytest.mark.integration
class TestAlignmentSubsampleIntegration:

    @patch("builtins.print")
    def test_basic_genes(self, mocked_print, tmp_path):
        gene_list = str(sample_dir / "gene_list_for_subsample.txt")
        prefix = str(tmp_path / "out")
        testargs = [
            "phykit",
            "alignment_subsample",
            "--mode", "genes",
            "-l", gene_list,
            "--number", "2",
            "--seed", "42",
            "-o", prefix,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        lines = _read_lines(f"{prefix}.txt")
        assert len(lines) == 2

    @patch("builtins.print")
    def test_basic_sites(self, mocked_print, tmp_path):
        aln = str(sample_dir / "test_alignment_0.fa")
        prefix = str(tmp_path / "out")
        testargs = [
            "phykit",
            "alignment_subsample",
            "--mode", "sites",
            "-a", aln,
            "--number", "5",
            "--seed", "42",
            "-o", prefix,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        seqs = _read_fasta(f"{prefix}.fa")
        assert len(seqs) == 4
        for seq in seqs.values():
            assert len(seq) == 5

    @patch("builtins.print")
    def test_basic_partitions(self, mocked_print, tmp_path):
        aln = str(sample_dir / "test_alignment_concat_123.fa")
        part = str(sample_dir / "test_alignment_concat_123.partition")
        prefix = str(tmp_path / "out")
        testargs = [
            "phykit",
            "alignment_subsample",
            "--mode", "partitions",
            "-a", aln,
            "-p", part,
            "--number", "2",
            "--seed", "42",
            "-o", prefix,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        seqs = _read_fasta(f"{prefix}.fa")
        # 2 partitions of 3 columns each = 6
        for seq in seqs.values():
            assert len(seq) == 6
        part_lines = _read_lines(f"{prefix}.partition")
        assert len(part_lines) == 2

    @patch("builtins.print")
    def test_alias_subsample(self, mocked_print, tmp_path):
        gene_list = str(sample_dir / "gene_list_for_subsample.txt")
        prefix = str(tmp_path / "out")
        testargs = [
            "phykit",
            "subsample",
            "--mode", "genes",
            "-l", gene_list,
            "--number", "2",
            "--seed", "42",
            "-o", prefix,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        lines = _read_lines(f"{prefix}.txt")
        assert len(lines) == 2

    @patch("builtins.print")
    def test_alias_aln_subsample(self, mocked_print, tmp_path):
        gene_list = str(sample_dir / "gene_list_for_subsample.txt")
        prefix = str(tmp_path / "out")
        testargs = [
            "phykit",
            "aln_subsample",
            "--mode", "genes",
            "-l", gene_list,
            "--number", "2",
            "--seed", "42",
            "-o", prefix,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        lines = _read_lines(f"{prefix}.txt")
        assert len(lines) == 2

    @patch("builtins.print")
    def test_json_output(self, mocked_print, tmp_path):
        gene_list = str(sample_dir / "gene_list_for_subsample.txt")
        prefix = str(tmp_path / "out")
        testargs = [
            "phykit",
            "alignment_subsample",
            "--mode", "genes",
            "-l", gene_list,
            "--number", "2",
            "--seed", "42",
            "-o", prefix,
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        # The JSON print call should be the only call (single line)
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["mode"] == "genes"
        assert payload["total"] == 3
        assert payload["selected"] == 2
        assert payload["bootstrap"] is False
        assert payload["seed"] == 42
        assert len(payload["output_files"]) == 1
