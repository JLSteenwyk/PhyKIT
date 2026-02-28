import json
import sys
from pathlib import Path

import pytest
from mock import patch

from phykit.phykit import Phykit


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
MULTI_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_multi_traits.tsv")


@pytest.mark.integration
class TestPhylogeneticPCAIntegration:
    @patch("builtins.print")
    def test_ppca_json(self, mocked_print):
        testargs = [
            "phykit",
            "phylogenetic_pca",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert "eigenvalues" in payload
        assert "proportion_of_variance" in payload
        assert "loadings" in payload
        assert "scores" in payload
        assert "PC1" in payload["eigenvalues"]

    @patch("builtins.print")
    def test_ppca_alias(self, mocked_print):
        testargs = [
            "phykit",
            "ppca",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Eigenvalues:" in all_output
        assert "Scores:" in all_output

    @patch("builtins.print")
    def test_ppca_lambda(self, mocked_print):
        testargs = [
            "phykit",
            "ppca",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "-m", "lambda",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert "lambda" in payload
        assert "log_likelihood" in payload
        assert "eigenvalues" in payload

    @patch("builtins.print")
    def test_ppca_corr_mode(self, mocked_print):
        testargs = [
            "phykit",
            "phyl_pca",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "--mode", "corr",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert "eigenvalues" in payload
        assert "scores" in payload

    @patch("builtins.print")
    def test_ppca_plot(self, mocked_print, tmp_path):
        plot_path = str(tmp_path / "pca_plot.png")
        testargs = [
            "phykit",
            "ppca",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "--plot",
            "--plot-output", plot_path,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert Path(plot_path).exists()
        assert Path(plot_path).stat().st_size > 0

    @patch("builtins.print")
    def test_ppca_plot_tree(self, mocked_print, tmp_path):
        plot_path = str(tmp_path / "pca_tree.png")
        testargs = [
            "phykit",
            "ppca",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "--plot",
            "--plot-tree",
            "--plot-output", plot_path,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert Path(plot_path).exists()
        assert Path(plot_path).stat().st_size > 0

    @patch("builtins.print")
    def test_ppca_color_by_column(self, mocked_print, tmp_path):
        plot_path = str(tmp_path / "pca_color.png")
        testargs = [
            "phykit",
            "ppca",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "--plot",
            "--color-by", "body_mass",
            "--plot-output", plot_path,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert Path(plot_path).exists()
        assert Path(plot_path).stat().st_size > 0

    @patch("builtins.print")
    def test_ppca_plot_tree_with_color_by(self, mocked_print, tmp_path):
        plot_path = str(tmp_path / "pca_both.png")
        testargs = [
            "phykit",
            "ppca",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "--plot",
            "--plot-tree",
            "--color-by", "body_mass",
            "--plot-output", plot_path,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert Path(plot_path).exists()
        assert Path(plot_path).stat().st_size > 0
