"""
Integration tests for phylo_heatmap (phylogenetic heatmap).

Analogous to R's phytools::phylo.heatmap(). Tests end-to-end CLI
execution with sample tree and multi-column trait data.
"""
from mock import patch
from pathlib import Path
import json
import sys

import pytest

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestPhyloHeatmap:
    @patch("builtins.print")
    def test_phylo_heatmap_creates_file(self, mocked_print, tmp_path):
        output = str(tmp_path / "heatmap.png")
        testargs = [
            "phykit",
            "phylo_heatmap",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/tree_simple_multi_traits.tsv",
            "-o", output,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    @patch("builtins.print")
    def test_phylo_heatmap_alias_pheatmap(self, mocked_print, tmp_path):
        output = str(tmp_path / "heatmap.png")
        testargs = [
            "phykit",
            "pheatmap",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/tree_simple_multi_traits.tsv",
            "-o", output,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    @patch("builtins.print")
    def test_phylo_heatmap_alias_ph(self, mocked_print, tmp_path):
        output = str(tmp_path / "heatmap.png")
        testargs = [
            "phykit",
            "ph",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/tree_simple_multi_traits.tsv",
            "-o", output,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    @patch("builtins.print")
    def test_phylo_heatmap_json(self, mocked_print, tmp_path):
        output = str(tmp_path / "heatmap.png")
        testargs = [
            "phykit",
            "phylo_heatmap",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/tree_simple_multi_traits.tsv",
            "-o", output,
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["n_taxa"] == 8
        assert payload["n_traits"] == 3

    @patch("builtins.print")
    def test_phylo_heatmap_with_standardize(self, mocked_print, tmp_path):
        output = str(tmp_path / "heatmap_std.png")
        testargs = [
            "phykit",
            "phylo_heatmap",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/tree_simple_multi_traits.tsv",
            "-o", output,
            "--standardize",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    @patch("builtins.print")
    def test_phylo_heatmap_custom_split(self, mocked_print, tmp_path):
        output = str(tmp_path / "heatmap_split.png")
        testargs = [
            "phykit",
            "phylo_heatmap",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/tree_simple_multi_traits.tsv",
            "-o", output,
            "--split", "0.5",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    @patch("builtins.print")
    def test_phylo_heatmap_pdf_output(self, mocked_print, tmp_path):
        output = str(tmp_path / "heatmap.pdf")
        testargs = [
            "phykit",
            "phylo_heatmap",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/tree_simple_multi_traits.tsv",
            "-o", output,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    @patch("builtins.print")
    def test_phylo_heatmap_circular(self, mocked_print, tmp_path):
        output = str(tmp_path / "heatmap_circular.png")
        testargs = [
            "phykit",
            "phylo_heatmap",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/tree_simple_multi_traits.tsv",
            "-o", output,
            "--circular",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()
        assert Path(output).stat().st_size > 0
