import json
import os
import tempfile

import pytest
import numpy as np
from argparse import Namespace
from pathlib import Path

from phykit.services.tree.phylo_impute import PhyloImpute
from phykit.errors import PhykitUserError


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
MISSING_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_traits_missing.tsv")
MULTI_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_multi_traits.tsv")


def _make_args(output_path, **overrides):
    defaults = dict(
        tree=TREE_SIMPLE,
        trait_data=MISSING_TRAITS_FILE,
        output=output_path,
        gene_trees=None,
        json=False,
    )
    defaults.update(overrides)
    return Namespace(**defaults)


def _build_service(output_path, **overrides):
    args = _make_args(output_path, **overrides)
    return PhyloImpute(args)


class TestPhyloImpute:
    def test_parse_na_values(self, tmp_path):
        """NA, ?, and empty values are all parsed as NaN."""
        out = str(tmp_path / "imputed.tsv")
        svc = _build_service(out)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits, missing_info = svc._parse_trait_file_with_na(
            MISSING_TRAITS_FILE, tree_tips
        )
        # bear has NA for brain_size, sea_lion has ? for longevity,
        # monkey has NA for longevity
        missing_taxa_traits = [(m["taxon"], m["trait"]) for m in missing_info]
        assert ("bear", "brain_size") in missing_taxa_traits
        assert ("sea_lion", "longevity") in missing_taxa_traits
        assert ("monkey", "longevity") in missing_taxa_traits

        # Check that the values are NaN
        assert np.isnan(traits["bear"][1])  # brain_size index 1
        assert np.isnan(traits["sea_lion"][2])  # longevity index 2
        assert np.isnan(traits["monkey"][2])  # longevity index 2

    def test_imputed_values_finite(self, tmp_path):
        """All imputed values should be finite numbers."""
        out = str(tmp_path / "imputed.tsv")
        svc = _build_service(out)
        svc.run()

        # Read output and verify no NaN/inf
        with open(out) as f:
            lines = f.readlines()
        for line in lines[1:]:
            parts = line.strip().split("\t")
            for val_str in parts[1:]:
                val = float(val_str)
                assert np.isfinite(val), f"Non-finite value {val_str} in output"

    def test_se_positive(self, tmp_path, capsys):
        """All standard errors should be positive."""
        out = str(tmp_path / "imputed.tsv")
        args = _make_args(out, json=True)
        svc = PhyloImpute(args)
        svc.run()

        captured = capsys.readouterr()
        payload = json.loads(captured.out)
        for entry in payload["imputed"]:
            assert entry["se"] > 0, (
                f"SE for {entry['taxon']}:{entry['trait']} = {entry['se']} <= 0"
            )

    def test_ci_contains_estimate(self, tmp_path, capsys):
        """CI lower < estimate < CI upper for all imputed values."""
        out = str(tmp_path / "imputed.tsv")
        args = _make_args(out, json=True)
        svc = PhyloImpute(args)
        svc.run()

        captured = capsys.readouterr()
        payload = json.loads(captured.out)
        for entry in payload["imputed"]:
            assert entry["ci_lower"] < entry["value"] < entry["ci_upper"], (
                f"CI [{entry['ci_lower']}, {entry['ci_upper']}] does not "
                f"contain estimate {entry['value']} for "
                f"{entry['taxon']}:{entry['trait']}"
            )

    def test_output_file_created(self, tmp_path):
        """Output TSV file should be created."""
        out = str(tmp_path / "imputed.tsv")
        svc = _build_service(out)
        svc.run()
        assert os.path.exists(out)

    def test_output_no_missing(self, tmp_path):
        """Output TSV should have no NA/missing values."""
        out = str(tmp_path / "imputed.tsv")
        svc = _build_service(out)
        svc.run()

        with open(out) as f:
            lines = f.readlines()

        missing_markers = {"NA", "na", "Na", "?", ""}
        for line_idx, line in enumerate(lines[1:], 2):
            parts = line.strip().split("\t")
            for i, val_str in enumerate(parts[1:]):
                assert val_str.strip() not in missing_markers, (
                    f"Missing value on line {line_idx}, column {i + 1}"
                )

    def test_complete_data_unchanged(self, tmp_path):
        """Taxa with no missing values should have unchanged values."""
        out = str(tmp_path / "imputed.tsv")
        svc = _build_service(out)
        svc.run()

        # Read original
        with open(MISSING_TRAITS_FILE) as f:
            orig_lines = f.readlines()
        orig_data = {}
        for line in orig_lines[1:]:
            parts = line.strip().split("\t")
            taxon = parts[0]
            vals = parts[1:]
            # Only track taxa with no missing values
            if all(v.strip() not in {"NA", "na", "Na", "?", ""} for v in vals):
                orig_data[taxon] = [float(v) for v in vals]

        # Read output
        with open(out) as f:
            out_lines = f.readlines()
        out_data = {}
        for line in out_lines[1:]:
            parts = line.strip().split("\t")
            out_data[parts[0]] = [float(v) for v in parts[1:]]

        for taxon, orig_vals in orig_data.items():
            for j, (orig_v, out_v) in enumerate(zip(orig_vals, out_data[taxon])):
                assert abs(orig_v - out_v) < 1e-4, (
                    f"Value changed for {taxon} trait {j}: "
                    f"{orig_v} -> {out_v}"
                )

    def test_json_output(self, tmp_path, capsys):
        """JSON output should have correct structure."""
        out = str(tmp_path / "imputed.tsv")
        args = _make_args(out, json=True)
        svc = PhyloImpute(args)
        svc.run()

        captured = capsys.readouterr()
        payload = json.loads(captured.out)

        assert "n_taxa" in payload
        assert "n_traits" in payload
        assert "trait_names" in payload
        assert "n_missing" in payload
        assert "vcv_type" in payload
        assert "imputed" in payload
        assert "output_file" in payload
        assert payload["n_traits"] == 3
        assert payload["n_missing"] == 3
        assert payload["vcv_type"] == "BM"
        assert len(payload["imputed"]) == 3

        for entry in payload["imputed"]:
            assert "taxon" in entry
            assert "trait" in entry
            assert "value" in entry
            assert "se" in entry
            assert "ci_lower" in entry
            assert "ci_upper" in entry

    def test_single_trait_works(self, tmp_path):
        """Single-trait imputation should work."""
        # Create a single-trait file with a missing value
        single_trait = tmp_path / "single_trait.tsv"
        single_trait.write_text(
            "taxon\tbody_mass\n"
            "raccoon\t1.04\n"
            "bear\tNA\n"
            "sea_lion\t2.30\n"
            "seal\t1.88\n"
            "monkey\t0.60\n"
            "cat\t0.56\n"
            "weasel\t-0.30\n"
            "dog\t1.18\n"
        )

        out = str(tmp_path / "imputed_single.tsv")
        args = _make_args(out, trait_data=str(single_trait))
        svc = PhyloImpute(args)
        svc.run()

        assert os.path.exists(out)
        with open(out) as f:
            lines = f.readlines()
        # Should have 9 lines (header + 8 taxa)
        assert len(lines) == 9

        # All values should be finite
        for line in lines[1:]:
            parts = line.strip().split("\t")
            val = float(parts[1])
            assert np.isfinite(val)

    def test_all_missing_one_trait(self, tmp_path):
        """If all taxa are missing one trait, use phylogenetic mean."""
        all_missing = tmp_path / "all_missing.tsv"
        all_missing.write_text(
            "taxon\tbody_mass\tnew_trait\n"
            "raccoon\t1.04\tNA\n"
            "bear\t2.39\tNA\n"
            "sea_lion\t2.30\tNA\n"
            "seal\t1.88\tNA\n"
            "monkey\t0.60\tNA\n"
            "cat\t0.56\tNA\n"
            "weasel\t-0.30\tNA\n"
            "dog\t1.18\tNA\n"
        )

        out = str(tmp_path / "imputed_all_missing.tsv")
        args = _make_args(out, trait_data=str(all_missing))
        svc = PhyloImpute(args)

        # This should fail because no taxa have complete data for all traits
        with pytest.raises(PhykitUserError):
            svc.run()

    def test_text_output_format(self, tmp_path, capsys):
        """Text output should contain expected fields."""
        out = str(tmp_path / "imputed.tsv")
        args = _make_args(out, json=False)
        svc = PhyloImpute(args)
        svc.run()

        captured = capsys.readouterr()
        text = captured.out
        assert "Phylogenetic Imputation" in text
        assert "Taxa:" in text
        assert "Traits:" in text
        assert "Missing values:" in text
        assert "VCV:" in text
        assert "Imputed values:" in text
        assert "Output:" in text
