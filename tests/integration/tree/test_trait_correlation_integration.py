import json
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
MULTI_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_multi_traits.tsv")


@pytest.mark.integration
class TestTraitCorrelationIntegration:
    def test_basic_invocation(self, tmp_path):
        out = str(tmp_path / "corr.png")
        testargs = [
            "phykit",
            "trait_correlation",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "-o", out,
        ]
        with patch("builtins.print") as mocked_print:
            with patch.object(sys, "argv", testargs):
                Phykit()
        assert os.path.exists(out)

    def test_alias_trait_corr(self, tmp_path):
        out = str(tmp_path / "corr.png")
        testargs = [
            "phykit",
            "trait_corr",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "-o", out,
        ]
        with patch("builtins.print") as mocked_print:
            with patch.object(sys, "argv", testargs):
                Phykit()
        assert os.path.exists(out)

    def test_alias_phylo_corr(self, tmp_path):
        out = str(tmp_path / "corr.png")
        testargs = [
            "phykit",
            "phylo_corr",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "-o", out,
        ]
        with patch("builtins.print") as mocked_print:
            with patch.object(sys, "argv", testargs):
                Phykit()
        assert os.path.exists(out)

    def test_custom_alpha(self, tmp_path):
        out = str(tmp_path / "corr.png")
        testargs = [
            "phykit",
            "trait_correlation",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "-o", out,
            "--alpha", "0.01",
        ]
        with patch("builtins.print") as mocked_print:
            with patch.object(sys, "argv", testargs):
                Phykit()
        assert os.path.exists(out)

    def test_cluster(self, tmp_path):
        out = str(tmp_path / "corr_cluster.png")
        testargs = [
            "phykit",
            "trait_correlation",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "-o", out,
            "--cluster",
        ]
        with patch("builtins.print") as mocked_print:
            with patch.object(sys, "argv", testargs):
                Phykit()
        assert os.path.exists(out)

    def test_json_output(self, tmp_path):
        out = str(tmp_path / "corr.png")
        testargs = [
            "phykit",
            "trait_correlation",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_FILE,
            "-o", out,
            "--json",
        ]
        with patch("builtins.print") as mocked_print:
            with patch.object(sys, "argv", testargs):
                Phykit()

        # Find the JSON call (the last print call)
        payload = json.loads(mocked_print.call_args.args[0])
        assert "n_taxa" in payload
        assert "n_traits" in payload
        assert "trait_names" in payload
        assert "correlation_matrix" in payload
        assert "p_value_matrix" in payload
        assert "significant_pairs" in payload
        assert payload["n_traits"] == 3
        assert payload["n_taxa"] == 8
        assert len(payload["trait_names"]) == 3
