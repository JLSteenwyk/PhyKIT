import json
import os
import sys
from pathlib import Path

import pytest
from mock import patch

from phykit.phykit import Phykit


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
MISSING_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_traits_missing.tsv")


@pytest.mark.integration
class TestPhyloImputeIntegration:
    def test_basic_invocation(self, tmp_path):
        out = str(tmp_path / "imputed.tsv")
        testargs = [
            "phykit",
            "phylo_impute",
            "-t", TREE_SIMPLE,
            "-d", MISSING_TRAITS_FILE,
            "-o", out,
        ]
        with patch("builtins.print") as mocked_print:
            with patch.object(sys, "argv", testargs):
                Phykit()
        assert os.path.exists(out)

    def test_alias_impute(self, tmp_path):
        out = str(tmp_path / "imputed.tsv")
        testargs = [
            "phykit",
            "impute",
            "-t", TREE_SIMPLE,
            "-d", MISSING_TRAITS_FILE,
            "-o", out,
        ]
        with patch("builtins.print") as mocked_print:
            with patch.object(sys, "argv", testargs):
                Phykit()
        assert os.path.exists(out)

    def test_alias_phylo_imp(self, tmp_path):
        out = str(tmp_path / "imputed.tsv")
        testargs = [
            "phykit",
            "phylo_imp",
            "-t", TREE_SIMPLE,
            "-d", MISSING_TRAITS_FILE,
            "-o", out,
        ]
        with patch("builtins.print") as mocked_print:
            with patch.object(sys, "argv", testargs):
                Phykit()
        assert os.path.exists(out)

    def test_json_output(self, tmp_path):
        out = str(tmp_path / "imputed.tsv")
        testargs = [
            "phykit",
            "phylo_impute",
            "-t", TREE_SIMPLE,
            "-d", MISSING_TRAITS_FILE,
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
        assert "n_missing" in payload
        assert "imputed" in payload
        assert payload["n_traits"] == 3
        assert payload["n_taxa"] == 8
        assert payload["n_missing"] == 3
        assert len(payload["imputed"]) == 3

    def test_output_file_complete(self, tmp_path):
        """Output TSV should be readable with no NAs."""
        out = str(tmp_path / "imputed.tsv")
        testargs = [
            "phykit",
            "phylo_impute",
            "-t", TREE_SIMPLE,
            "-d", MISSING_TRAITS_FILE,
            "-o", out,
        ]
        with patch("builtins.print"):
            with patch.object(sys, "argv", testargs):
                Phykit()

        with open(out) as f:
            lines = f.readlines()

        # Header + 8 taxa
        assert len(lines) == 9

        missing_markers = {"NA", "na", "Na", "?", ""}
        for line in lines[1:]:
            parts = line.strip().split("\t")
            assert len(parts) == 4  # taxon + 3 traits
            for val_str in parts[1:]:
                assert val_str.strip() not in missing_markers
                val = float(val_str)
                assert not (val != val)  # not NaN
