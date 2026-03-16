"""
Integration tests for quartet_pie (quartet concordance pie chart).

Ground truth gCF/gDF values computed via manual bipartition matching
on tests/sample_files/tree_simple.tre + gene_trees_simple.nwk.
Values cross-checked against concordance_asr output.
"""
from mock import patch
from pathlib import Path
import json
import sys

import pytest

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestQuartetPie:
    @patch("builtins.print")
    def test_quartet_pie_creates_file(self, mocked_print, tmp_path):
        output = str(tmp_path / "qpie.png")
        testargs = [
            "phykit",
            "quartet_pie",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-g", f"{here.parent.parent.parent}/sample_files/gene_trees_simple.nwk",
            "-o", output,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    @patch("builtins.print")
    def test_quartet_pie_json_ground_truth(self, mocked_print, tmp_path):
        output = str(tmp_path / "qpie.png")
        testargs = [
            "phykit",
            "quartet_pie",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-g", f"{here.parent.parent.parent}/sample_files/gene_trees_simple.nwk",
            "-o", output,
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["n_taxa"] == 8
        assert payload["n_gene_trees"] == 10
        assert payload["input_mode"] == "native"

        # Verify ground truth: bear+raccoon fully concordant
        br = next(n for n in payload["nodes"]
                  if set(n["node_tips"]) == {"bear", "raccoon"})
        assert br["gCF"] == 1.0
        assert br["concordant_count"] == 7

        # monkey+cat+weasel: gCF ~ 0.47, gDF1 ~ 0.53
        mcw = next(n for n in payload["nodes"]
                   if set(n["node_tips"]) == {"monkey", "cat", "weasel"})
        assert mcw["gCF"] == pytest.approx(9 / 19, abs=0.01)
        assert mcw["gDF1"] == pytest.approx(10 / 19, abs=0.01)

    @patch("builtins.print")
    def test_quartet_pie_alias_qpie(self, mocked_print, tmp_path):
        output = str(tmp_path / "qpie.png")
        testargs = [
            "phykit",
            "qpie",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-g", f"{here.parent.parent.parent}/sample_files/gene_trees_simple.nwk",
            "-o", output,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    @patch("builtins.print")
    def test_quartet_pie_pdf_output(self, mocked_print, tmp_path):
        output = str(tmp_path / "qpie.pdf")
        testargs = [
            "phykit",
            "quartet_pie",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-g", f"{here.parent.parent.parent}/sample_files/gene_trees_simple.nwk",
            "-o", output,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    @patch("builtins.print")
    def test_quartet_pie_with_annotate(self, mocked_print, tmp_path):
        output = str(tmp_path / "qpie_annot.png")
        testargs = [
            "phykit",
            "quartet_pie",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-g", f"{here.parent.parent.parent}/sample_files/gene_trees_simple.nwk",
            "-o", output,
            "--annotate",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    @patch("builtins.print")
    def test_quartet_pie_wastral_support3(self, mocked_print, tmp_path):
        """wASTRAL --support 3 output is parsed in ASTRAL mode."""
        output = str(tmp_path / "qpie_wastral.png")
        testargs = [
            "phykit",
            "quartet_pie",
            "-t", f"{here.parent.parent.parent}/sample_files/wastral_support3.tre",
            "-o", output,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    @patch("builtins.print")
    def test_quartet_pie_wastral_json(self, mocked_print, tmp_path):
        """wASTRAL --support 3 JSON output has correct q-values."""
        output = str(tmp_path / "qpie_wastral.png")
        testargs = [
            "phykit",
            "quartet_pie",
            "-t", f"{here.parent.parent.parent}/sample_files/wastral_support3.tre",
            "-o", output,
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["input_mode"] == "astral"
        assert payload["n_taxa"] == 7
        assert len(payload["nodes"]) == 6

        # Verify the node with most discordance (Vallesia + Haplophyton)
        vh = next(
            n for n in payload["nodes"]
            if len(n["node_tips"]) == 2
            and any("Vallesia" in t for t in n["node_tips"])
        )
        assert vh["gCF"] == pytest.approx(0.6449, abs=0.01)
        assert vh["gDF1"] == pytest.approx(0.2190, abs=0.01)
        assert vh["gDF2"] == pytest.approx(0.1361, abs=0.01)
