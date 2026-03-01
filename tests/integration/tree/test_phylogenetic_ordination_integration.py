import json
import sys
from pathlib import Path

import pytest
from mock import patch

# Import umap before any print patching to avoid numba/mock conflict
import umap  # noqa: F401

from phykit.phykit import Phykit


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
MULTI_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_multi_traits.tsv")


@pytest.mark.integration
class TestPhylogeneticOrdinationIntegration:
    # --- PCA tests (backwards compat) ---

    @patch("builtins.print")
    def test_pca_json(self, mocked_print):
        testargs = [
            "phykit",
            "phylogenetic_ordination",
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
            "--correction", "lambda",
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

    # --- Dimreduce tests (backwards compat) ---

    @patch("builtins.print")
    def test_tsne_default(self, mocked_print):
        testargs = [
            "phykit",
            "phylogenetic_ordination",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "--method", "tsne",
            "--seed", "42",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Method: tsne" in all_output
        assert "Embedding:" in all_output

    @patch("builtins.print")
    def test_umap_invocation(self, mocked_print):
        testargs = [
            "phykit",
            "phylogenetic_ordination",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "--method", "umap",
            "--seed", "42",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Method: umap" in all_output
        assert "Embedding:" in all_output

    @patch("builtins.print")
    def test_alias_dimreduce(self, mocked_print):
        testargs = [
            "phykit",
            "dimreduce",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "--method", "tsne",
            "--seed", "42",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Embedding:" in all_output

    @patch("builtins.print")
    def test_alias_pdr(self, mocked_print):
        testargs = [
            "phykit",
            "pdr",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "--method", "tsne",
            "--seed", "42",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Embedding:" in all_output

    @patch("builtins.print")
    def test_dimreduce_json_output(self, mocked_print):
        testargs = [
            "phykit",
            "dimreduce",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "--method", "tsne",
            "--seed", "42",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["method"] == "tsne"
        assert "embedding" in payload
        assert "parameters" in payload

    @patch("builtins.print")
    def test_dimreduce_lambda_correction(self, mocked_print):
        testargs = [
            "phykit",
            "dimreduce",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "--method", "tsne",
            "--correction", "lambda",
            "--seed", "42",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert "lambda" in payload
        assert "log_likelihood" in payload

    @patch("builtins.print")
    def test_custom_perplexity(self, mocked_print):
        testargs = [
            "phykit",
            "dimreduce",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "--method", "tsne",
            "--perplexity", "2.0",
            "--seed", "42",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["parameters"]["perplexity"] == 2.0

    @patch("builtins.print")
    def test_dimreduce_plot_creation(self, mocked_print, tmp_path):
        plot_path = str(tmp_path / "dimreduce_plot.png")
        testargs = [
            "phykit",
            "dimreduce",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "--method", "tsne",
            "--seed", "42",
            "--plot",
            "--plot-output", plot_path,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert Path(plot_path).exists()
        assert Path(plot_path).stat().st_size > 0

    # --- New ordination aliases ---

    @patch("builtins.print")
    def test_alias_ordination(self, mocked_print):
        testargs = [
            "phykit",
            "ordination",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert "eigenvalues" in payload

    @patch("builtins.print")
    def test_alias_ord(self, mocked_print):
        testargs = [
            "phykit",
            "ord",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "--method", "umap",
            "--seed", "42",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["method"] == "umap"
        assert "embedding" in payload
