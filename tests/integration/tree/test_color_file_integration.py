"""
Integration tests for --color-file support across all 11 phylogram-drawing services.

Each test creates a small color annotation file in tmp_path, runs the service
with --color-file, and asserts the output plot file was created.
"""
import os
import sys
import tempfile
from pathlib import Path

import pytest
from mock import patch

from phykit.phykit import Phykit

here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
TREE_CHARMAP = str(SAMPLE_FILES / "tree_character_map.tre")
TREE_OTHER = str(SAMPLE_FILES / "tree_simple_other_topology.tre")
GENE_TREES = str(SAMPLE_FILES / "gene_trees_simple.nwk")
TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_traits.tsv")
MULTI_TRAITS = str(SAMPLE_FILES / "tree_simple_multi_traits.tsv")
DISCRETE_TRAITS = str(SAMPLE_FILES / "tree_simple_discrete_traits.tsv")
CHARMATRIX = str(SAMPLE_FILES / "character_matrix_simple.tsv")
REGIMES_FILE = str(SAMPLE_FILES / "tree_simple_regimes.tsv")


def _write_color_file_simple(tmp_path):
    """Color file compatible with tree_simple.tre (8 taxa)."""
    p = tmp_path / "colors.tsv"
    p.write_text(
        "bear\tlabel\t#ff0000\n"
        "bear,raccoon\trange\t#ffe0e0\tUrsids\n"
    )
    return str(p)


def _write_color_file_charmap(tmp_path):
    """Color file compatible with tree_character_map.tre (6 taxa: A-F)."""
    p = tmp_path / "colors.tsv"
    p.write_text(
        "A\tlabel\t#ff0000\n"
        "A,B\trange\t#ffe0e0\tGroupAB\n"
    )
    return str(p)


@pytest.mark.integration
class TestColorFileIntegration:
    # 1. quartet_pie
    @patch("builtins.print")
    def test_quartet_pie_color_file(self, mocked_print, tmp_path):
        color_file = _write_color_file_simple(tmp_path)
        output = str(tmp_path / "qpie_colored.png")
        testargs = [
            "phykit", "quartet_pie",
            "-t", TREE_SIMPLE,
            "-g", GENE_TREES,
            "-o", output,
            "--color-file", color_file,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    # 2. character_map
    @patch("builtins.print")
    def test_character_map_color_file(self, mocked_print, tmp_path):
        color_file = _write_color_file_charmap(tmp_path)
        output = str(tmp_path / "charmap_colored.png")
        testargs = [
            "phykit", "character_map",
            "-t", TREE_CHARMAP,
            "-d", CHARMATRIX,
            "-o", output,
            "--color-file", color_file,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    # 3. phylo_heatmap
    @patch("builtins.print")
    def test_phylo_heatmap_color_file(self, mocked_print, tmp_path):
        color_file = _write_color_file_simple(tmp_path)
        output = str(tmp_path / "heatmap_colored.png")
        testargs = [
            "phykit", "phylo_heatmap",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS,
            "-o", output,
            "--color-file", color_file,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    # 4. cophylo
    @patch("builtins.print")
    def test_cophylo_color_file(self, mocked_print, tmp_path):
        color_file = _write_color_file_simple(tmp_path)
        output = str(tmp_path / "cophylo_colored.png")
        testargs = [
            "phykit", "cophylo",
            "-t", TREE_SIMPLE,
            "-t2", TREE_OTHER,
            "-o", output,
            "--color-file", color_file,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    # 5. ancestral_reconstruction (continuous contMap)
    @patch("builtins.print")
    def test_ancestral_reconstruction_contmap_color_file(self, mocked_print, tmp_path):
        color_file = _write_color_file_simple(tmp_path)
        output = str(tmp_path / "asr_contmap_colored.png")
        testargs = [
            "phykit", "asr",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
            "--plot", output,
            "--color-file", color_file,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    # 5b. ancestral_reconstruction (discrete)
    @patch("builtins.print")
    def test_ancestral_reconstruction_discrete_color_file(self, mocked_print, tmp_path):
        color_file = _write_color_file_simple(tmp_path)
        output = str(tmp_path / "asr_discrete_colored.png")
        testargs = [
            "phykit", "asr",
            "-t", TREE_SIMPLE,
            "-d", DISCRETE_TRAITS,
            "--type", "discrete",
            "--plot", output,
            "--color-file", color_file,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    # 6. cont_map
    @patch("builtins.print")
    def test_cont_map_color_file(self, mocked_print, tmp_path):
        color_file = _write_color_file_simple(tmp_path)
        output = str(tmp_path / "contmap_colored.png")
        testargs = [
            "phykit", "cont_map",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
            "-o", output,
            "--color-file", color_file,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    # 7. density_map
    @patch("builtins.print")
    def test_density_map_color_file(self, mocked_print, tmp_path):
        color_file = _write_color_file_simple(tmp_path)
        output = str(tmp_path / "densitymap_colored.png")
        testargs = [
            "phykit", "density_map",
            "-t", TREE_SIMPLE,
            "-d", DISCRETE_TRAITS,
            "-c", "diet",
            "-n", "5",
            "--seed", "42",
            "-o", output,
            "--color-file", color_file,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    # 8. stochastic_character_map
    @patch("builtins.print")
    def test_stochastic_character_map_color_file(self, mocked_print, tmp_path):
        color_file = _write_color_file_simple(tmp_path)
        output = str(tmp_path / "scm_colored.png")
        testargs = [
            "phykit", "stochastic_character_map",
            "-t", TREE_SIMPLE,
            "-d", DISCRETE_TRAITS,
            "-c", "diet",
            "-n", "5",
            "--seed", "42",
            "--plot", output,
            "--color-file", color_file,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    # 9. discordance_asymmetry
    @patch("builtins.print")
    def test_discordance_asymmetry_color_file(self, mocked_print, tmp_path):
        color_file = _write_color_file_simple(tmp_path)
        output = str(tmp_path / "dasym_colored.png")
        testargs = [
            "phykit", "discordance_asymmetry",
            "-t", TREE_SIMPLE,
            "-g", GENE_TREES,
            "--plot", output,
            "--color-file", color_file,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    # 10. rate_heterogeneity
    @patch("builtins.print")
    def test_rate_heterogeneity_color_file(self, mocked_print, tmp_path):
        color_file = _write_color_file_simple(tmp_path)
        output = str(tmp_path / "ratehet_colored.png")
        testargs = [
            "phykit", "rate_heterogeneity",
            "-t", TREE_SIMPLE,
            "-d", TRAITS_FILE,
            "-r", REGIMES_FILE,
            "--plot", output,
            "--color-file", color_file,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    # 11. concordance_asr
    @patch("builtins.print")
    def test_concordance_asr_color_file(self, mocked_print, tmp_path):
        color_file = _write_color_file_simple(tmp_path)
        output = str(tmp_path / "conc_asr_colored.png")
        testargs = [
            "phykit", "concordance_asr",
            "-t", TREE_SIMPLE,
            "-g", GENE_TREES,
            "-d", TRAITS_FILE,
            "--plot", output,
            "--color-file", color_file,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()
