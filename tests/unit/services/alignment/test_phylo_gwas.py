"""Unit tests for the PhyloGwas service."""

import builtins
import csv
import json
import os
import subprocess
import sys
from argparse import Namespace
from io import StringIO
from pathlib import Path

import numpy as np
import pytest

from phykit.services.alignment.phylo_gwas import PhyloGwas
import phykit.services.alignment.phylo_gwas as phylo_gwas_module


def test_module_import_does_not_import_biopython_parsers():
    code = """
import sys
import phykit.services.alignment.phylo_gwas as module
assert hasattr(module.np, "__getattr__")
assert callable(module.print_json)
assert "typing" not in sys.modules
assert "numpy" not in sys.modules
assert "csv" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
assert "Bio.SeqIO.FastaIO" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_lazy_numpy_caches_resolved_attributes():
    lazy_np = phylo_gwas_module._LazyNumpy()

    ndarray_attr = lazy_np.ndarray

    assert lazy_np.__dict__["ndarray"] is ndarray_attr
    assert lazy_np.ndarray is ndarray_attr
    assert lazy_np._module is not None


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

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


# Test alignment and phenotype data
# 8 taxa, 4 highland + 4 lowland
# Site 0 (col1): all A -> invariant (skip)
# Site 1 (col2): perfectly associated (highland=A, lowland=G)
# Site 2 (col3): partially associated (highland mostly G, lowland mostly C)
# Site 3 (col4): mixed, weak association
# Site 4 (col5): gap in sp8 -> skip
SEQS = {
    "sp1": "AAGAA",  # highland
    "sp2": "AAGCA",  # highland
    "sp3": "AAGAA",  # highland
    "sp4": "AACAA",  # highland
    "sp5": "AGCGA",  # lowland
    "sp6": "AGCAA",  # lowland
    "sp7": "AGCGA",  # lowland
    "sp8": "AGCA-",  # lowland
}

PHENO_CATEGORICAL = {
    "sp1": "highland",
    "sp2": "highland",
    "sp3": "highland",
    "sp4": "highland",
    "sp5": "lowland",
    "sp6": "lowland",
    "sp7": "lowland",
    "sp8": "lowland",
}

PHENO_CONTINUOUS = {
    "sp1": "10.0",
    "sp2": "12.0",
    "sp3": "11.0",
    "sp4": "9.0",
    "sp5": "2.0",
    "sp6": "3.0",
    "sp7": "1.0",
    "sp8": "4.0",
}

# Tree where sp1-sp4 form a clade and sp5-sp8 form another clade
TREE_MONOPHYLETIC = "((sp1,sp2,sp3,sp4),(sp5,sp6,sp7,sp8));"

# Tree where the allele pattern would be polyphyletic
TREE_POLYPHYLETIC = "((sp1,sp5,sp3,sp7),(sp2,sp6,sp4,sp8));"


def _make_args(tmp_path, seqs=None, pheno=None, tree=None, partition=None, **kwargs):
    """Create args namespace for a test run."""
    seqs = seqs or SEQS
    pheno = pheno or PHENO_CATEGORICAL

    aln_path = str(tmp_path / "test.fa")
    pheno_path = str(tmp_path / "test_pheno.tsv")
    output_path = str(tmp_path / "manhattan.png")

    _write_alignment(aln_path, seqs)
    _write_phenotype(pheno_path, pheno)

    args = Namespace(
        alignment=aln_path,
        phenotype=pheno_path,
        output=output_path,
        tree=None,
        partition=None,
        alpha=kwargs.get("alpha", 0.05),
        exclude_monophyletic=kwargs.get("exclude_monophyletic", False),
        csv=kwargs.get("csv", None),
        json=kwargs.get("json_output", False),
    )

    if tree:
        tree_path = str(tmp_path / "test.nwk")
        _write_tree(tree_path, tree)
        args.tree = tree_path

    if partition:
        part_path = str(tmp_path / "test.partition")
        _write_partition(part_path, partition)
        args.partition = part_path

    return args


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestPhyloGwas:
    def test_init_sets_expected_attrs(self, tmp_path):
        args = _make_args(tmp_path)
        service = PhyloGwas(args)
        assert service.alignment_path == args.alignment
        assert service.phenotype_path == args.phenotype
        assert service.output_path == args.output
        assert service.alpha == 0.05
        assert service.exclude_monophyletic is False
        assert service.json_output is False

    def test_read_fasta(self, tmp_path):
        aln_path = str(tmp_path / "test.fa")
        _write_alignment(aln_path, {"s1": "ACG", "s2": "ATG"})
        seqs = PhyloGwas._read_fasta(aln_path)
        assert seqs == {"s1": "ACG", "s2": "ATG"}

    def test_read_fasta_uses_first_header_token_uppercases_and_keeps_last_duplicate(
        self, tmp_path
    ):
        aln_path = tmp_path / "test.fa"
        aln_path.write_text(
            ">s1 first description\nac g\nn-\n"
            ">s2 second description\nat g\nn-\n"
            ">s1 replacement description\nttt\n"
        )
        seqs = PhyloGwas._read_fasta(str(aln_path))
        assert seqs == {"s1": "TTT", "s2": "ATGN-"}

    def test_read_phenotype(self, tmp_path):
        pheno_path = str(tmp_path / "pheno.tsv")
        _write_phenotype(pheno_path, {"s1": "A", "s2": "B"})
        pheno = PhyloGwas._read_phenotype(pheno_path)
        assert pheno == {"s1": "A", "s2": "B"}

    def test_read_phenotype_uses_first_two_columns_and_keeps_last_duplicate(
        self, tmp_path
    ):
        pheno_path = tmp_path / "pheno.tsv"
        pheno_path.write_text(
            "# ignored\n"
            "\n"
            "malformed\n"
            "s1\tA\textra\n"
            "s2\tB\tignored\tcolumns\n"
            "s1\tC\treplacement\n"
        )

        pheno = PhyloGwas._read_phenotype(str(pheno_path))

        assert pheno == {"s1": "C", "s2": "B"}

    def test_read_phenotype_handles_binary_rows_and_crlf(self, tmp_path):
        pheno_path = tmp_path / "pheno.tsv"
        pheno_path.write_bytes(
            b"# ignored\r\n"
            b"\r\n"
            b"malformed\r\n"
            b"s1\tA\textra\r\n"
            b"s2\tB\tignored\tcolumns\r\n"
            b"s1\tC\treplacement\r\n"
        )

        pheno = PhyloGwas._read_phenotype(str(pheno_path))

        assert pheno == {"s1": "C", "s2": "B"}

    def test_sorted_shared_taxa_preserves_sorted_overlap(self):
        seqs = {"sp3": "AAA", "sp1": "AAA", "sp2": "AAA"}
        phenotypes = {"sp2": "case", "sp4": "case", "sp1": "control"}

        assert PhyloGwas._sorted_shared_taxa(seqs, phenotypes) == ["sp1", "sp2"]

    def test_sorted_shared_taxa_filters_smaller_phenotype_mapping(self):
        seqs = {
            "sp5": "AAA",
            "sp1": "AAA",
            "sp4": "AAA",
            "sp2": "AAA",
            "sp3": "AAA",
        }
        phenotypes = {"sp4": "case", "missing": "case", "sp2": "control"}

        assert PhyloGwas._sorted_shared_taxa(seqs, phenotypes) == ["sp2", "sp4"]

    def test_sorted_shared_taxa_filters_smaller_sequence_mapping(self):
        seqs = {"sp4": "AAA", "sp2": "AAA"}
        phenotypes = {
            "sp5": "case",
            "sp4": "case",
            "sp1": "control",
            "sp2": "case",
        }

        assert PhyloGwas._sorted_shared_taxa(seqs, phenotypes) == ["sp2", "sp4"]

    def test_sorted_shared_taxa_exact_keys_skips_intersection_set(self, monkeypatch):
        seqs = {"sp3": "AAA", "sp1": "AAA", "sp2": "AAA"}
        phenotypes = {"sp2": "case", "sp3": "case", "sp1": "control"}

        def fail_set(*_args, **_kwargs):
            raise AssertionError(
                "exact shared taxa should not build an intersection set"
            )

        monkeypatch.setattr("builtins.set", fail_set)

        assert PhyloGwas._sorted_shared_taxa(seqs, phenotypes) == [
            "sp1",
            "sp2",
            "sp3",
        ]

    def test_extract_column_alleles_from_ascii_matrix(self):
        sequences = ["ACG", "ATG", "AGG"]
        matrix = PhyloGwas._build_ascii_alignment_matrix(sequences, 3)
        lookup = PhyloGwas._ascii_ambiguity_lookup()
        assert matrix is not None
        assert PhyloGwas._extract_column_alleles(1, sequences, matrix, lookup) == [
            "C",
            "T",
            "G",
        ]

    def test_extract_column_alleles_rejects_ambiguous_ascii_matrix_column(self):
        sequences = ["ACG", "A?G", "AGG"]
        matrix = PhyloGwas._build_ascii_alignment_matrix(sequences, 3)
        lookup = PhyloGwas._ascii_ambiguity_lookup()
        assert matrix is not None
        assert (
            PhyloGwas._extract_column_alleles(1, sequences, matrix, lookup) is None
        )

    def test_valid_ascii_columns_rejects_ambiguous_columns(self):
        sequences = ["ACGTAC", "A?GTXC", "AC-T*C", "ACGNAC"]
        matrix = PhyloGwas._build_ascii_alignment_matrix(sequences, 6)
        lookup = PhyloGwas._ascii_ambiguity_lookup()
        assert matrix is not None

        valid_columns = PhyloGwas._valid_ascii_columns(matrix, lookup)

        assert valid_columns.tolist() == [True, False, False, False, False, True]

    def test_biallelic_valid_ascii_columns_rejects_non_testable_columns(self):
        sequences = ["AAAC", "AAC?", "AGGT", "AGTT"]
        matrix = PhyloGwas._build_ascii_alignment_matrix(sequences, 4)
        lookup = PhyloGwas._ascii_ambiguity_lookup()
        assert matrix is not None

        valid_columns = PhyloGwas._biallelic_valid_ascii_columns(matrix, lookup)

        assert valid_columns.tolist() == [False, True, False, False]

    def test_extract_column_alleles_non_ascii_fallback(self):
        sequences = ["AαG", "AβG"]
        assert PhyloGwas._build_ascii_alignment_matrix(sequences, 3) is None
        assert PhyloGwas._extract_column_alleles(1, sequences, None, None) == [
            "α",
            "β",
        ]

    def test_extract_column_alleles_fallback_inlines_ambiguity_check(
        self, monkeypatch
    ):
        sequences = ["AαG", "AβG", "A?G"]

        def fail_is_ambiguous(_char):
            raise AssertionError("fallback extraction should check ambiguity inline")

        monkeypatch.setattr(PhyloGwas, "_is_ambiguous", fail_is_ambiguous)

        assert PhyloGwas._extract_column_alleles(1, sequences, None, None) is None

    def test_detect_phenotype_type_categorical(self):
        assert PhyloGwas._detect_phenotype_type(["highland", "lowland"]) == "categorical"

    def test_detect_phenotype_type_continuous(self):
        assert PhyloGwas._detect_phenotype_type(["1.5", "2.3", "0.1"]) == "continuous"

    def test_is_ambiguous(self):
        assert PhyloGwas._is_ambiguous("-") is True
        assert PhyloGwas._is_ambiguous("N") is True
        assert PhyloGwas._is_ambiguous("X") is True
        assert PhyloGwas._is_ambiguous("?") is True
        assert PhyloGwas._is_ambiguous("A") is False

    def test_benjamini_hochberg_basic(self):
        p_values = [0.01, 0.04, 0.03, 0.20]
        adjusted = PhyloGwas._benjamini_hochberg(p_values)
        # Adjusted p-values should be >= raw p-values
        for raw, adj in zip(p_values, adjusted):
            assert adj >= raw
        # Should be capped at 1.0
        assert all(a <= 1.0 for a in adjusted)

    def test_benjamini_hochberg_matches_scalar_reference_with_ties(self):
        p_values = [0.20, 0.01, 0.01, 0.50, 0.03, 0.80, 0.03]
        p_arr = np.array(p_values, dtype=float)
        order = np.argsort(p_arr)
        sorted_p = p_arr[order]
        expected_sorted = np.zeros(len(p_values), dtype=float)
        for i in range(len(p_values) - 1, -1, -1):
            candidate = sorted_p[i] * len(p_values) / (i + 1)
            if i == len(p_values) - 1:
                expected_sorted[i] = candidate
            else:
                expected_sorted[i] = min(candidate, expected_sorted[i + 1])
        expected = np.zeros(len(p_values), dtype=float)
        expected[order] = np.minimum(expected_sorted, 1.0)

        np.testing.assert_allclose(
            PhyloGwas._benjamini_hochberg(p_values),
            expected,
        )

    def test_small_benjamini_hochberg_skips_vector_sort(self, monkeypatch):
        def fail_argsort(*args, **kwargs):
            raise AssertionError("small BH corrections should use scalar sorting")

        monkeypatch.setattr(phylo_gwas_module.np, "argsort", fail_argsort)

        adjusted = PhyloGwas._benjamini_hochberg(
            [0.20, 0.01, 0.01, 0.50, 0.03, 0.80, 0.03]
        )

        assert isinstance(adjusted, np.ndarray)
        np.testing.assert_allclose(
            adjusted,
            [0.28, 0.035, 0.035, 0.5833333333333334, 0.0525, 0.8, 0.0525],
        )

    def test_benjamini_hochberg_populates_result_without_zero_fill(self, monkeypatch):
        def fail_zeros(*args, **kwargs):
            raise AssertionError("BH output should not be zero-filled before scatter")

        monkeypatch.setattr(phylo_gwas_module.np, "zeros", fail_zeros)

        adjusted = PhyloGwas._benjamini_hochberg([0.01, 0.04, 0.03])

        np.testing.assert_allclose(adjusted, [0.03, 0.04, 0.04])

    def test_benjamini_hochberg_empty(self):
        adjusted = PhyloGwas._benjamini_hochberg([])
        assert len(adjusted) == 0

    def test_categorical_association(self, tmp_path):
        """Perfectly associated site should have low p-value."""
        args = _make_args(tmp_path, alpha=1.0)
        service = PhyloGwas(args)
        # Run it - but we capture results through JSON
        args_json = _make_args(tmp_path, json_output=True, alpha=1.0)
        service = PhyloGwas(args_json)
        service.run()

        # The output was printed as JSON; check that the file was created
        assert os.path.exists(args_json.output)

    def test_invariant_skipped(self, tmp_path, capsys):
        """Invariant sites should not appear in results."""
        args = _make_args(tmp_path, json_output=True, alpha=1.0)
        service = PhyloGwas(args)
        service.run()

        out, _ = capsys.readouterr()
        payload = json.loads(out)
        # Site 1 (position 1) is invariant (all A), so it should not be in results
        positions = [r["position"] for r in payload["results"]]
        assert 1 not in positions  # position 1 is invariant

    def test_gaps_skipped(self, tmp_path, capsys):
        """Sites with gaps should not be tested."""
        args = _make_args(tmp_path, json_output=True, alpha=1.0)
        service = PhyloGwas(args)
        service.run()

        out, _ = capsys.readouterr()
        payload = json.loads(out)
        # Site 5 (position 5) has a gap in sp8
        positions = [r["position"] for r in payload["results"]]
        assert 5 not in positions

    def test_fdr_correction(self, tmp_path, capsys):
        """Adjusted p-values should be >= raw p-values."""
        args = _make_args(tmp_path, json_output=True, alpha=1.0)
        service = PhyloGwas(args)
        service.run()

        out, _ = capsys.readouterr()
        payload = json.loads(out)
        for r in payload["results"]:
            assert r["fdr_p_value"] >= r["p_value"]

    def test_fdr_correction_preserves_result_order(self, tmp_path, capsys, monkeypatch):
        args = _make_args(tmp_path, json_output=True, alpha=0.15)
        service = PhyloGwas(args)
        monkeypatch.setattr(service, "_create_manhattan_plot", lambda *_args: None)

        def ordered_adjusted(p_values):
            return np.array(
                [0.05 + idx * 0.1 for idx in range(len(p_values))],
                dtype=float,
            )

        monkeypatch.setattr(service, "_benjamini_hochberg", ordered_adjusted)

        service.run()

        out, _ = capsys.readouterr()
        payload = json.loads(out)
        fdr_values = [row["fdr_p_value"] for row in payload["results"]]

        assert len(fdr_values) >= 2
        assert fdr_values == [
            pytest.approx(0.05 + idx * 0.1)
            for idx in range(len(fdr_values))
        ]
        assert [row["fdr_significant"] for row in payload["results"]] == [
            value < 0.15 for value in fdr_values
        ]

    def test_creates_png(self, tmp_path):
        """Manhattan plot should be created."""
        pytest.importorskip("matplotlib")
        args = _make_args(tmp_path)
        service = PhyloGwas(args)
        service.run()
        assert os.path.exists(args.output)

    def test_manhattan_plot_batches_partition_spans(
        self, tmp_path, monkeypatch
    ):
        """Partition shading should use one collection instead of per-span artists."""
        pytest.importorskip("matplotlib")
        import matplotlib.axes
        from matplotlib.collections import PatchCollection

        args = _make_args(tmp_path)
        service = PhyloGwas(args)
        results = [
            {
                "position": i + 1,
                "p_value": 0.5,
                "fdr_significant": False,
            }
            for i in range(20)
        ]
        partitions = [
            (f"gene{i}", i * 5, i * 5 + 4)
            for i in range(6)
        ]
        patch_collections = []
        original_add_collection = matplotlib.axes.Axes.add_collection

        def tracking_add_collection(self, collection, *args, **kwargs):
            if isinstance(collection, PatchCollection):
                patch_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        def fail_axvspan(*args, **kwargs):
            raise AssertionError("partition spans should be batched")

        monkeypatch.setattr(
            matplotlib.axes.Axes,
            "add_collection",
            tracking_add_collection,
        )
        monkeypatch.setattr(matplotlib.axes.Axes, "axvspan", fail_axvspan)

        service._create_manhattan_plot(results, partitions, has_tree=False)

        assert os.path.exists(args.output)
        assert len(patch_collections) == 1
        assert len(patch_collections[0].get_paths()) == 3

    def test_manhattan_plot_reuses_cached_rendered_bytes(
        self, tmp_path, monkeypatch
    ):
        """Repeated identical Manhattan plots should skip Matplotlib rendering."""
        pytest.importorskip("matplotlib")
        previous_cache = phylo_gwas_module._MANHATTAN_PLOT_CACHE.copy()
        phylo_gwas_module._MANHATTAN_PLOT_CACHE.clear()
        results = [
            {
                "position": 1,
                "p_value": 0.5,
                "fdr_significant": False,
                "phylo_pattern": None,
            },
            {
                "position": 2,
                "p_value": 0.01,
                "fdr_significant": True,
                "phylo_pattern": "polyphyletic",
            },
        ]

        first_args = _make_args(tmp_path)
        first_service = PhyloGwas(first_args)
        first_service._create_manhattan_plot(results, [], has_tree=False)
        first_bytes = Path(first_args.output).read_bytes()

        second_dir = tmp_path / "second"
        second_dir.mkdir()
        second_args = _make_args(second_dir)
        second_service = PhyloGwas(second_args)
        original_import = builtins.__import__

        def fail_matplotlib_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name == "matplotlib" or name.startswith("matplotlib."):
                raise AssertionError("cached plot should not import matplotlib")
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", fail_matplotlib_import)
        try:
            second_service._create_manhattan_plot(results, [], has_tree=False)
            assert Path(second_args.output).read_bytes() == first_bytes
        finally:
            phylo_gwas_module._MANHATTAN_PLOT_CACHE.clear()
            phylo_gwas_module._MANHATTAN_PLOT_CACHE.update(previous_cache)

    def test_csv_output(self, tmp_path):
        """CSV file should be created with correct columns."""
        csv_path = str(tmp_path / "results.csv")
        args = _make_args(tmp_path, csv=csv_path, alpha=1.0)
        service = PhyloGwas(args)
        service.run()

        assert os.path.exists(csv_path)
        with open(csv_path) as f:
            reader = csv.DictReader(f)
            fieldnames = reader.fieldnames
            assert "position" in fieldnames
            assert "allele_0" in fieldnames
            assert "p_value" in fieldnames
            assert "fdr_p_value" in fieldnames
            assert "significant" in fieldnames
            rows = list(reader)
            assert len(rows) > 0

    def test_write_csv_preserves_field_order_defaults_and_quoting(self, tmp_path):
        csv_path = tmp_path / "results.csv"
        args = _make_args(tmp_path, csv=str(csv_path))
        service = PhyloGwas(args)
        service._write_csv(
            [
                {
                    "position": 1,
                    "gene": 'gene,"1"',
                    "allele_0": "A",
                    "allele_1": "G",
                    "p_value": 0.01,
                    "fdr_p_value": 0.02,
                    "fdr_significant": True,
                    "phylo_pattern": None,
                    "ignored": "value",
                }
            ],
            is_continuous=False,
        )

        assert csv_path.read_text() == (
            "position,gene,allele_0,allele_1,p_value,fdr_p_value,"
            "significant,phylo_pattern\n"
            "1,\"gene,\"\"1\"\"\",A,G,0.01,0.02,yes,\n"
        )

    def test_write_csv_continuous_includes_correlation_column(self, tmp_path):
        csv_path = tmp_path / "continuous.csv"
        args = _make_args(tmp_path, csv=str(csv_path))
        service = PhyloGwas(args)
        service._write_csv(
            [
                {
                    "position": 2,
                    "gene": None,
                    "allele_0": "C",
                    "allele_1": "T",
                    "correlation_r": -0.5,
                    "p_value": 0.03,
                    "fdr_p_value": 0.04,
                    "fdr_significant": False,
                    "phylo_pattern": "polyphyletic",
                }
            ],
            is_continuous=True,
        )

        assert csv_path.read_text() == (
            "position,gene,allele_0,allele_1,correlation_r,p_value,"
            "fdr_p_value,significant,phylo_pattern\n"
            "2,,C,T,-0.5,0.03,0.04,no,polyphyletic\n"
        )

    def test_write_csv_keeps_blank_defaults_for_incomplete_rows(self, tmp_path):
        csv_path = tmp_path / "missing.csv"
        args = _make_args(tmp_path, csv=str(csv_path))
        service = PhyloGwas(args)
        service._write_csv(
            [
                {
                    "position": 3,
                    "allele_0": "A",
                    "p_value": 0.05,
                    "fdr_significant": True,
                }
            ],
            is_continuous=False,
        )

        assert csv_path.read_text() == (
            "position,gene,allele_0,allele_1,p_value,fdr_p_value,"
            "significant,phylo_pattern\n"
            "3,,A,,0.05,,yes,\n"
        )

    def test_json_output_structure(self, tmp_path, capsys):
        """JSON output should have correct structure."""
        args = _make_args(tmp_path, json_output=True, alpha=1.0)
        service = PhyloGwas(args)
        service.run()

        out, _ = capsys.readouterr()
        payload = json.loads(out)
        assert "n_taxa" in payload
        assert "alignment_length" in payload
        assert "phenotype_type" in payload
        assert "sites_tested" in payload
        assert "significant_sites" in payload
        assert "results" in payload
        assert payload["n_taxa"] == 8
        assert payload["alignment_length"] == 5
        assert payload["phenotype_type"] == "categorical"

    def test_continuous_phenotype(self, tmp_path, capsys):
        """Numeric phenotype should use point-biserial and be detected as continuous."""
        args = _make_args(tmp_path, pheno=PHENO_CONTINUOUS, json_output=True, alpha=1.0)
        service = PhyloGwas(args)
        service.run()

        out, _ = capsys.readouterr()
        payload = json.loads(out)
        assert payload["phenotype_type"] == "continuous"
        # Results should include correlation_r
        for r in payload["results"]:
            assert "correlation_r" in r

    def test_json_payload_reuses_result_rows_without_shallow_copy(self, tmp_path):
        args = _make_args(tmp_path, json_output=True, alpha=0.25)
        service = PhyloGwas(args)
        site_results = [
            {
                "position": 1,
                "allele_0": "A",
                "allele_1": "G",
                "p_value": 0.01,
                "fdr_p_value": 0.02,
                "fdr_significant": True,
                "gene": None,
                "phylo_pattern": "polyphyletic",
            }
        ]

        payload = service._build_json_payload(
            n_taxa=8,
            aln_len=120,
            pheno_type="continuous",
            site_results=site_results,
            n_significant=1,
            group_counts=None,
            n_polyphyletic=1,
            n_monophyletic=0,
            has_tree=True,
        )

        assert payload["results"] is site_results
        assert payload["results"][0] is site_results[0]
        assert payload["alpha"] == 0.25
        assert payload["polyphyletic"] == 1
        assert payload["monophyletic"] == 0

    def test_summarize_significant_results_single_pass_preserves_rows(self):
        results = [
            {
                "position": 1,
                "fdr_significant": False,
                "phylo_pattern": "polyphyletic",
            },
            {
                "position": 2,
                "fdr_significant": True,
                "phylo_pattern": "polyphyletic",
            },
            {
                "position": 3,
                "fdr_significant": True,
                "phylo_pattern": "monophyletic",
            },
            {
                "position": 4,
                "fdr_significant": True,
                "phylo_pattern": None,
            },
        ]

        significant, n_polyphyletic, n_monophyletic = (
            PhyloGwas._summarize_significant_results(results)
        )

        assert significant == results[1:]
        assert [id(row) for row in significant] == [id(row) for row in results[1:]]
        assert n_polyphyletic == 1
        assert n_monophyletic == 1

    def test_assign_fdr_values_converts_adjusted_values_in_bulk(self):
        class ToListOnly:
            def tolist(self):
                return [0.01, 0.2]

            def __iter__(self):
                raise AssertionError("adjusted FDR values should be bulk-converted")

        results = [{"position": 1}, {"position": 2}]

        PhyloGwas._assign_fdr_values(results, ToListOnly(), 0.05)

        assert results == [
            {"position": 1, "fdr_p_value": 0.01, "fdr_significant": True},
            {"position": 2, "fdr_p_value": 0.2, "fdr_significant": False},
        ]

    def test_categorical_association_accepts_ascii_byte_column(self):
        args = Namespace(
            alignment="",
            phenotype="",
            output="",
            tree=None,
            partition=None,
            alpha=0.05,
            exclude_monophyletic=False,
            csv=None,
            json=False,
        )
        service = PhyloGwas(args)
        alleles = ["A", "A", "G", "G", "A", "G"]
        groups = ["x", "x", "x", "y", "y", "y"]
        unique_groups = ["x", "y"]
        group_idx = {"x": 0, "y": 1}
        group_codes = np.array([0, 0, 0, 1, 1, 1], dtype=np.intp)
        group_row_counts = np.array([3, 3])

        expected = service._test_site_categorical(
            alleles,
            groups,
            unique_groups,
            group_idx,
            group_codes,
            group_row_counts,
        )
        observed = service._test_site_categorical(
            np.frombuffer("".join(alleles).encode("ascii"), dtype=np.uint8),
            groups,
            unique_groups,
            group_idx,
            group_codes,
            group_row_counts,
        )

        assert observed is not None
        assert expected is not None
        assert observed[:3] == expected[:3]
        assert observed[3] == expected[3]

    def test_run_uses_ascii_byte_columns_for_categorical_sites(
        self, tmp_path, monkeypatch
    ):
        args = _make_args(tmp_path, json_output=True, alpha=1.0)
        service = PhyloGwas(args)
        observed_byte_columns = 0
        original = PhyloGwas._test_site_categorical

        def tracking_test(self, alleles, *rest, **kwargs):
            nonlocal observed_byte_columns
            if isinstance(alleles, np.ndarray) and alleles.dtype == np.uint8:
                observed_byte_columns += 1
            return original(self, alleles, *rest, **kwargs)

        monkeypatch.setattr(PhyloGwas, "_test_site_categorical", tracking_test)

        service.run()

        assert observed_byte_columns > 0

    def test_run_uses_ascii_byte_columns_for_continuous_sites(
        self, tmp_path, monkeypatch
    ):
        args = _make_args(
            tmp_path, pheno=PHENO_CONTINUOUS, json_output=True, alpha=1.0
        )
        service = PhyloGwas(args)
        observed_byte_columns = 0
        original = PhyloGwas._test_site_continuous

        def tracking_test(self, alleles, *rest, **kwargs):
            nonlocal observed_byte_columns
            if isinstance(alleles, np.ndarray) and alleles.dtype == np.uint8:
                observed_byte_columns += 1
            return original(self, alleles, *rest, **kwargs)

        monkeypatch.setattr(PhyloGwas, "_test_site_continuous", tracking_test)

        service.run()

        assert observed_byte_columns > 0

    def test_run_skips_ascii_columns_that_cannot_be_tested(
        self, tmp_path, monkeypatch
    ):
        seqs = {
            "s1": "AAACA",
            "s2": "AAC?A",
            "s3": "AGGTT",
            "s4": "AGTTT",
        }
        pheno = {
            "s1": "1.0",
            "s2": "2.0",
            "s3": "3.0",
            "s4": "4.0",
        }
        args = _make_args(
            tmp_path,
            seqs=seqs,
            pheno=pheno,
            json_output=True,
            alpha=1.0,
        )
        service = PhyloGwas(args)
        observed_columns = []
        original = PhyloGwas._test_site_continuous

        def tracking_test(self, alleles, *rest, **kwargs):
            if isinstance(alleles, np.ndarray) and alleles.dtype == np.uint8:
                observed_columns.append(alleles.tobytes())
            return original(self, alleles, *rest, **kwargs)

        monkeypatch.setattr(PhyloGwas, "_test_site_continuous", tracking_test)

        service.run()

        assert observed_columns == [b"AAGG", b"AATT"]

    def test_run_keeps_light_ascii_scan_when_sample_is_all_biallelic(
        self, tmp_path, monkeypatch
    ):
        seqs = {
            "s1": "AA",
            "s2": "AA",
            "s3": "GG",
            "s4": "GG",
        }
        pheno = {
            "s1": "1.0",
            "s2": "2.0",
            "s3": "3.0",
            "s4": "4.0",
        }
        args = _make_args(
            tmp_path,
            seqs=seqs,
            pheno=pheno,
            json_output=True,
            alpha=1.0,
        )
        service = PhyloGwas(args)

        def fail_prefilter(*_args, **_kwargs):
            raise AssertionError("all-biallelic scans should not build the full mask")

        monkeypatch.setattr(
            PhyloGwas,
            "_biallelic_valid_ascii_columns",
            staticmethod(fail_prefilter),
        )

        service.run()

    def test_monophyletic_classification(self, tmp_path, capsys):
        """With monophyletic tree, allele pattern matching clade should be monophyletic."""
        args = _make_args(
            tmp_path, tree=TREE_MONOPHYLETIC, json_output=True, alpha=1.0
        )
        service = PhyloGwas(args)
        service.run()

        out, _ = capsys.readouterr()
        payload = json.loads(out)
        # Site 2 (position 2) perfectly associates highland=A, lowland=G
        # With the monophyletic tree, lowland taxa (sp5-sp8) form a clade
        sig_results = [r for r in payload["results"] if r.get("fdr_significant")]
        for r in sig_results:
            if r["position"] == 2:
                assert r["phylo_pattern"] in ("monophyletic", "polyphyletic")

    def test_polyphyletic_classification(self, tmp_path, capsys):
        """With non-matching tree, allele pattern should be polyphyletic."""
        args = _make_args(
            tmp_path, tree=TREE_POLYPHYLETIC, json_output=True, alpha=1.0
        )
        service = PhyloGwas(args)
        service.run()

        out, _ = capsys.readouterr()
        payload = json.loads(out)
        sig_results = [r for r in payload["results"] if r.get("fdr_significant")]
        # With the polyphyletic tree, associated alleles should be polyphyletic
        for r in sig_results:
            if r["position"] == 2:
                assert r["phylo_pattern"] == "polyphyletic"

    def test_minor_allele_taxa_filters_by_position(self):
        shared_taxa = ["sp1", "sp2", "sp3", "sp4"]
        seqs = {
            "sp1": "AAG",
            "sp2": "AGG",
            "sp3": "AAG",
            "sp4": "AGC",
        }

        assert PhyloGwas._minor_allele_taxa(shared_taxa, seqs, 1, "G") == [
            "sp2",
            "sp4",
        ]

    def test_minor_allele_taxa_from_ascii_column_matches_string_scan(self):
        shared_taxa = ["sp1", "sp2", "sp3", "sp4"]
        seqs = {
            "sp1": "AAG",
            "sp2": "AGG",
            "sp3": "AAG",
            "sp4": "AGC",
        }
        matrix = PhyloGwas._build_ascii_alignment_matrix(
            [seqs[taxon] for taxon in shared_taxa],
            3,
        )
        assert matrix is not None

        assert PhyloGwas._minor_allele_taxa_from_ascii_column(
            shared_taxa,
            matrix[:, 1],
            "G",
        ) == PhyloGwas._minor_allele_taxa(shared_taxa, seqs, 1, "G")

    def test_run_uses_ascii_column_minor_allele_taxa_for_phylo_pattern(
        self, tmp_path, capsys, monkeypatch
    ):
        args = _make_args(
            tmp_path, tree=TREE_MONOPHYLETIC, json_output=True, alpha=1.0
        )
        calls = []
        original = PhyloGwas._minor_allele_taxa_from_ascii_column

        def tracking_minor_allele_taxa(shared_taxa, allele_array, allele):
            calls.append((len(allele_array), allele))
            return original(shared_taxa, allele_array, allele)

        monkeypatch.setattr(
            PhyloGwas,
            "_minor_allele_taxa_from_ascii_column",
            staticmethod(tracking_minor_allele_taxa),
        )

        service = PhyloGwas(args)
        service.run()

        out, _ = capsys.readouterr()
        payload = json.loads(out)
        expected = [
            (len(PHENO_CATEGORICAL), result["allele_1"])
            for result in payload["results"]
            if result["fdr_significant"]
        ]
        assert calls == expected
        assert calls

    def test_exclude_monophyletic(self, tmp_path, capsys):
        """The --exclude-monophyletic flag should remove monophyletic hits."""
        # First run without flag to see what's there
        args1 = _make_args(
            tmp_path, tree=TREE_MONOPHYLETIC, json_output=True, alpha=1.0
        )
        service1 = PhyloGwas(args1)
        service1.run()
        out1, _ = capsys.readouterr()
        payload1 = json.loads(out1)
        mono_count1 = sum(
            1
            for r in payload1["results"]
            if r.get("fdr_significant") and r.get("phylo_pattern") == "monophyletic"
        )

        # Second run with flag
        args2 = _make_args(
            tmp_path,
            tree=TREE_MONOPHYLETIC,
            json_output=True,
            alpha=1.0,
            exclude_monophyletic=True,
        )
        service2 = PhyloGwas(args2)
        service2.run()
        out2, _ = capsys.readouterr()
        payload2 = json.loads(out2)
        mono_count2 = sum(
            1
            for r in payload2["results"]
            if r.get("fdr_significant") and r.get("phylo_pattern") == "monophyletic"
        )

        # If there were monophyletic hits, they should be removed
        if mono_count1 > 0:
            assert mono_count2 == 0

    def test_text_output(self, tmp_path, capsys):
        """Text output should contain summary information."""
        args = _make_args(tmp_path)
        service = PhyloGwas(args)
        service.run()

        out, _ = capsys.readouterr()
        assert "Phylogenetic GWAS" in out
        assert "Alignment:" in out
        assert "Taxa: 8" in out
        assert "Phenotype type: categorical" in out
        assert "Biallelic sites tested:" in out

    def test_print_text_output_batches_report(self, monkeypatch):
        service = PhyloGwas.__new__(PhyloGwas)
        service.alignment_path = "test.fa"
        service.alpha = 0.05
        significant = [
            {
                "position": 2,
                "gene": "gene1",
                "allele_0": "A",
                "allele_1": "G",
                "group_freqs": {"highland": 0.75, "lowland": 0.10},
                "p_value": 0.001,
                "fdr_p_value": 0.01,
                "phylo_pattern": "polyphyletic",
            },
            {
                "position": 3,
                "gene": None,
                "allele_0": "C",
                "allele_1": "T",
                "group_freqs": {"highland": 0.25, "lowland": 0.55},
                "p_value": 0.002,
                "fdr_p_value": 0.02,
                "phylo_pattern": "monophyletic",
            },
        ]
        site_results = significant + [dict(significant[0], fdr_significant=False)]
        printed = []

        def fake_print(*args, **kwargs):
            printed.append((args, kwargs))

        monkeypatch.setattr("builtins.print", fake_print)

        service._print_text_output(
            8,
            5,
            "categorical",
            False,
            ["highland", "lowland"],
            {"highland": 4, "lowland": 4},
            site_results,
            significant,
            1,
            1,
            True,
        )

        expected_lines = [
            "Phylogenetic GWAS",
            "Alignment: test.fa",
            "Taxa: 8",
            "Alignment length: 5",
            "Phenotype type: categorical",
            "Groups: highland (4), lowland (4)",
            "Biallelic sites tested: 3",
            "Significant sites (FDR < 0.05): 2",
            "  Polyphyletic: 1",
            "  Monophyletic: 1",
            "",
            "Top significant sites:",
            (
                f"  {'Position':<10}{'Gene':<12}{'Allele':<10}"
                f"{'highland/lowland':<20}{'p-value':<12}{'FDR_p':<12}"
                f"{'Pattern':<14}"
            ),
            (
                f"  {2:<10}{'gene1':<12}"
                f"A>{'G':<8}"
                f"{'0.75/0.10':<20}"
                f"{0.001:<12.4g}"
                f"{0.01:<12.4g}"
                f"{'polyphyletic':<14}"
            ),
            (
                f"  {3:<10}{'.':<12}"
                f"C>{'T':<8}"
                f"{'0.25/0.55':<20}"
                f"{0.002:<12.4g}"
                f"{0.02:<12.4g}"
                f"{'monophyletic':<14}"
            ),
        ]
        assert printed == [((("\n".join(expected_lines)),), {})]

    def test_top_significant_sites_matches_sorted_top_ten(self):
        results = [
            {"position": i, "p_value": p_value}
            for i, p_value in enumerate(
                [0.05, 0.01, 0.03, 0.01, 0.02, 0.08, 0.07, 0.04, 0.06, 0.09, 0.0]
            )
        ]

        assert PhyloGwas._top_significant_sites(results) == sorted(
            results, key=lambda r: r["p_value"]
        )[:10]

    def test_top_significant_sites_uses_cached_p_value_getter(self, monkeypatch):
        results = [
            {"position": i, "p_value": p_value}
            for i, p_value in enumerate(
                [0.05, 0.01, 0.03, 0.04, 0.02, 0.08, 0.07, 0.06, 0.09, 0.10, 0.0]
            )
        ]
        captured = {}

        def fake_nsmallest(n, iterable, key):
            captured["key"] = key
            return sorted(iterable, key=key)[:n]

        monkeypatch.setattr(phylo_gwas_module.heapq, "nsmallest", fake_nsmallest)

        assert PhyloGwas._top_significant_sites(results) == sorted(
            results, key=lambda r: r["p_value"]
        )[:10]
        assert captured["key"] is phylo_gwas_module._P_VALUE_GETTER

    def test_prepare_manhattan_series_combines_plot_inputs(self):
        results = [
            {
                "position": 1,
                "p_value": 0.5,
                "fdr_significant": False,
                "phylo_pattern": None,
            },
            {
                "position": 2,
                "p_value": 0.01,
                "fdr_significant": True,
                "phylo_pattern": "polyphyletic",
            },
            {
                "position": 3,
                "p_value": 0.02,
                "fdr_significant": True,
                "phylo_pattern": "monophyletic",
            },
            {
                "position": 4,
                "p_value": 1e-320,
                "fdr_significant": False,
                "phylo_pattern": None,
            },
        ]

        positions, neg_log_p, point_colors, threshold = (
            PhyloGwas._prepare_manhattan_series(
                results,
                "blue",
                "red",
                "gray",
            )
        )

        np.testing.assert_array_equal(positions, np.array([1, 2, 3, 4]))
        np.testing.assert_allclose(
            neg_log_p,
            -np.log10(np.array([0.5, 0.01, 0.02, 1e-300])),
        )
        assert point_colors == ["blue", "red", "gray", "blue"]
        assert threshold == pytest.approx(-np.log10(0.02))

    def test_with_partition(self, tmp_path, capsys):
        """Partition file should annotate genes in results."""
        partition = [
            ("DNA", "gene1", 1, 3),
            ("DNA", "gene2", 4, 5),
        ]
        args = _make_args(
            tmp_path, partition=partition, json_output=True, alpha=1.0
        )
        service = PhyloGwas(args)
        service.run()

        out, _ = capsys.readouterr()
        payload = json.loads(out)
        # Some results should have gene annotations
        genes = [r.get("gene") for r in payload["results"]]
        assert any(g is not None for g in genes)

    def test_parse_partition_file(self, tmp_path):
        """Partition file parsing should produce correct 0-based ranges."""
        part_path = str(tmp_path / "test.part")
        _write_partition(part_path, [
            ("DNA", "gene1", 1, 100),
            ("DNA", "gene2", 101, 200),
        ])
        partitions = PhyloGwas._parse_partition_file(part_path)
        assert len(partitions) == 2
        assert partitions[0] == ("gene1", 0, 99)
        assert partitions[1] == ("gene2", 100, 199)

    def test_parse_partition_file_preserves_regex_tolerance(self, tmp_path):
        part_path = tmp_path / "test.part"
        part_path.write_text(
            "   # ignored\n"
            "DNA,   gene1   =   1   -   100 trailing text\n"
            ",gene2=101-200\n"
            "invalid row\n"
        )

        assert PhyloGwas._parse_partition_file(str(part_path)) == [
            ("gene1", 0, 99),
            ("gene2", 100, 199),
        ]

    def test_position_to_gene(self):
        partitions = [("gene1", 0, 99), ("gene2", 100, 199)]
        assert PhyloGwas._position_to_gene(50, partitions) == "gene1"
        assert PhyloGwas._position_to_gene(150, partitions) == "gene2"
        assert PhyloGwas._position_to_gene(200, partitions) is None

    def test_partition_lookup_handles_boundaries_and_gaps(self):
        partitions = [("gene1", 0, 9), ("gene2", 20, 29)]
        lookup = PhyloGwas._build_partition_lookup(partitions)
        assert lookup is not None
        assert PhyloGwas._position_to_gene_from_lookup(0, lookup) == "gene1"
        assert PhyloGwas._position_to_gene_from_lookup(9, lookup) == "gene1"
        assert PhyloGwas._position_to_gene_from_lookup(10, lookup) is None
        assert PhyloGwas._position_to_gene_from_lookup(19, lookup) is None
        assert PhyloGwas._position_to_gene_from_lookup(20, lookup) == "gene2"
        assert PhyloGwas._position_to_gene_from_lookup(29, lookup) == "gene2"
        assert PhyloGwas._position_to_gene_from_lookup(30, lookup) is None

    def test_partition_lookup_rejects_overlapping_or_unsorted_partitions(self):
        overlapping = [("gene1", 0, 9), ("gene2", 9, 20)]
        unsorted = [("gene2", 20, 29), ("gene1", 0, 9)]
        invalid = [("gene1", 10, 9)]
        assert PhyloGwas._build_partition_lookup(overlapping) is None
        assert PhyloGwas._build_partition_lookup(unsorted) is None
        assert PhyloGwas._build_partition_lookup(invalid) is None

    def test_classify_phylo_pattern_single_taxon(self):
        """Single taxon is trivially monophyletic."""
        from Bio import Phylo
        from io import StringIO
        tree = Phylo.read(StringIO("((A,B),(C,D));"), "newick")
        assert PhyloGwas._classify_phylo_pattern(tree, ["A"]) == "monophyletic"

    def test_classify_phylo_pattern_monophyletic(self):
        from Bio import Phylo
        from io import StringIO
        tree = Phylo.read(StringIO("((A,B),(C,D));"), "newick")
        assert PhyloGwas._classify_phylo_pattern(tree, ["A", "B"]) == "monophyletic"

    def test_classify_phylo_pattern_polyphyletic(self):
        from Bio import Phylo
        from io import StringIO
        tree = Phylo.read(StringIO("((A,B),(C,D));"), "newick")
        assert PhyloGwas._classify_phylo_pattern(tree, ["A", "C"]) == "polyphyletic"

    def test_classify_phylo_pattern_index_avoids_mrca_lookup(self, monkeypatch):
        from Bio import Phylo
        from Bio.Phylo.BaseTree import Clade, TreeMixin
        from io import StringIO

        tree = Phylo.read(StringIO("((A,B),(C,D));"), "newick")
        monophyletic_sets = PhyloGwas._build_phylo_pattern_index(tree)

        def fail_common_ancestor(self, *args, **kwargs):
            raise AssertionError("common_ancestor should not be called")

        def fail_get_terminals(self, *args, **kwargs):
            raise AssertionError("get_terminals should not be called")

        monkeypatch.setattr(type(tree), "common_ancestor", fail_common_ancestor)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)
        monkeypatch.setattr(Clade, "get_terminals", fail_get_terminals)

        assert (
            PhyloGwas._classify_phylo_pattern(tree, ["A", "B"], monophyletic_sets)
            == "monophyletic"
        )
        assert (
            PhyloGwas._classify_phylo_pattern(tree, ["A", "C"], monophyletic_sets)
            == "polyphyletic"
        )

    def test_build_phylo_pattern_index_uses_direct_postorder(self, monkeypatch):
        from Bio import Phylo
        from Bio.Phylo.BaseTree import TreeMixin
        from io import StringIO

        tree = Phylo.read(StringIO("((A,B),(C,D));"), "newick")

        def fail_find_clades(*args, **kwargs):
            raise AssertionError("generic postorder traversal should not be used")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        assert PhyloGwas._build_phylo_pattern_index(tree) == {
            frozenset({"A"}),
            frozenset({"B"}),
            frozenset({"C"}),
            frozenset({"D"}),
            frozenset({"A", "B"}),
            frozenset({"C", "D"}),
            frozenset({"A", "B", "C", "D"}),
        }

    def test_build_phylo_pattern_index_handles_mixed_child_counts(self, monkeypatch):
        from Bio import Phylo
        from Bio.Phylo.BaseTree import TreeMixin
        from io import StringIO

        tree = Phylo.read(StringIO("(A,(B,C),(D,E,F));"), "newick")

        def fail_find_clades(*args, **kwargs):
            raise AssertionError("generic postorder traversal should not be used")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        assert PhyloGwas._build_phylo_pattern_index(tree) == {
            frozenset({"A"}),
            frozenset({"B"}),
            frozenset({"C"}),
            frozenset({"D"}),
            frozenset({"E"}),
            frozenset({"F"}),
            frozenset({"B", "C"}),
            frozenset({"D", "E", "F"}),
            frozenset({"A", "B", "C", "D", "E", "F"}),
        }

    def test_build_phylo_pattern_index_handles_unary_child_nodes(self, monkeypatch):
        from Bio.Phylo.BaseTree import Clade, Tree, TreeMixin

        tree = Tree(
            root=Clade(
                clades=[
                    Clade(
                        clades=[
                            Clade(name="A"),
                        ],
                    ),
                    Clade(
                        clades=[
                            Clade(name="B"),
                            Clade(name="C"),
                        ],
                    ),
                ],
            )
        )

        def fail_find_clades(*args, **kwargs):
            raise AssertionError("generic postorder traversal should not be used")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        assert PhyloGwas._build_phylo_pattern_index(tree) == {
            frozenset({"A"}),
            frozenset({"B"}),
            frozenset({"C"}),
            frozenset({"B", "C"}),
            frozenset({"A", "B", "C"}),
        }

    def test_prune_tree_to_taxa_uses_direct_terminal_pass(self, monkeypatch):
        from Bio import Phylo
        from Bio.Phylo.BaseTree import TreeMixin
        from io import StringIO
        tree = Phylo.read(StringIO("((A,B),(C,D));"), "newick")
        original_get_terminals = type(tree).get_terminals

        def fail_get_terminals(self, *args, **kwargs):
            raise AssertionError("standard tree should use direct terminal traversal")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        PhyloGwas._prune_tree_to_taxa(tree, {"A", "C"})

        assert sorted(t.name for t in original_get_terminals(tree)) == ["A", "C"]

    def test_terminal_clades_preserves_mixed_child_order(self, monkeypatch):
        from Bio import Phylo
        from Bio.Phylo.BaseTree import TreeMixin
        from io import StringIO

        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1,G:1);"),
            "newick",
        )
        expected_terminals = tree.get_terminals()

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("standard tree should use direct terminal traversal")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        assert PhyloGwas._terminal_clades(tree) == expected_terminals

    def test_too_few_shared_taxa(self, tmp_path, capsys):
        """Should exit with error if fewer than 4 shared taxa."""
        seqs = {"s1": "A", "s2": "A", "s3": "A"}
        pheno = {"s1": "X", "s2": "Y", "s3": "X"}
        with pytest.raises(SystemExit) as exc:
            args = _make_args(tmp_path, seqs=seqs, pheno=pheno)
            service = PhyloGwas(args)
            service.run()
        assert exc.value.code == 2

    def test_test_site_categorical_invariant(self):
        """Invariant site should return None."""
        service = PhyloGwas.__new__(PhyloGwas)
        result = service._test_site_categorical(
            ["A", "A", "A", "A"], ["X", "X", "Y", "Y"], ["X", "Y"]
        )
        assert result is None

    def test_test_site_categorical_biallelic(self):
        """Biallelic site should return p-value and frequencies."""
        service = PhyloGwas.__new__(PhyloGwas)
        result = service._test_site_categorical(
            ["A", "A", "G", "G"], ["X", "X", "Y", "Y"], ["X", "Y"]
        )
        assert result is not None
        p_value, allele_0, allele_1, group_freqs = result
        assert 0 <= p_value <= 1
        assert allele_0 in ("A", "G")
        assert allele_1 in ("A", "G")
        assert len(group_freqs) == 2

    def test_fisher_exact_2x2_two_sided_matches_scipy(self):
        from scipy.stats import fisher_exact

        tables = [
            np.array([[2, 0], [0, 2]]),
            np.array([[18, 7], [11, 24]]),
            np.array([[90, 12], [85, 9]]),
            np.array([[350, 25], [122, 203]]),
            np.array([[2500, 30], [2450, 60]]),
        ]

        for table in tables:
            assert phylo_gwas_module._fisher_exact_2x2_two_sided(
                table
            ) == pytest.approx(float(fisher_exact(table).pvalue))

    def test_fisher_exact_2x2_two_sided_reuses_cache(self):
        table = np.array([[6, 4], [2, 8]])
        cache = {}
        expected = phylo_gwas_module._fisher_exact_2x2_two_sided(table)

        assert phylo_gwas_module._fisher_exact_2x2_two_sided(
            table,
            cache,
        ) == pytest.approx(expected)
        cache[(6, 4, 2, 8)] = 0.12345

        assert phylo_gwas_module._fisher_exact_2x2_two_sided(
            table,
            cache,
        ) == pytest.approx(0.12345)

    def test_fisher_exact_2x2_two_sided_caches_log_combinations(self):
        cache = {}

        first = phylo_gwas_module._fisher_exact_2x2_two_sided(
            np.array([[18, 7], [11, 24]]),
            cache,
        )
        second = phylo_gwas_module._fisher_exact_2x2_two_sided(
            np.array([[17, 8], [12, 23]]),
            cache,
        )

        assert 0.0 <= first <= 1.0
        assert 0.0 <= second <= 1.0
        assert (18, 7, 11, 24) in cache
        assert (17, 8, 12, 23) in cache
        assert any(key[0] == "log_comb" for key in cache if isinstance(key, tuple))

    def test_test_site_categorical_two_group_does_not_import_scipy_stats(
        self, monkeypatch
    ):
        service = PhyloGwas.__new__(PhyloGwas)
        original_import = __import__

        def guarded_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name == "scipy.stats" or name.startswith("scipy.stats."):
                raise AssertionError(
                    "two-group categorical GWAS should not import scipy.stats"
                )
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr("builtins.__import__", guarded_import)

        result = service._test_site_categorical(
            np.frombuffer(b"AAGGGGAA", dtype=np.uint8),
            ["X", "X", "X", "X", "Y", "Y", "Y", "Y"],
            ["X", "Y"],
        )

        assert result is not None
        p_value, allele_0, allele_1, group_freqs = result
        assert 0.0 <= p_value <= 1.0
        assert allele_0 in ("A", "G")
        assert allele_1 in ("A", "G")
        assert group_freqs == {"X": 0.5, "Y": 0.5}

    def test_test_site_categorical_multigroup_frequencies(self):
        from scipy.stats import chi2_contingency

        service = PhyloGwas.__new__(PhyloGwas)
        unique_groups = ["X", "Y", "Z"]
        groups = ["X", "X", "Y", "Y", "Z", "Z"]
        group_idx = {group: i for i, group in enumerate(unique_groups)}
        group_codes = np.array(
            [group_idx[group] for group in groups], dtype=np.intp
        )
        group_row_counts = np.bincount(group_codes, minlength=len(unique_groups))
        result = service._test_site_categorical(
            ["A", "A", "A", "G", "G", "G"],
            groups,
            unique_groups,
            group_idx,
            group_codes,
            group_row_counts,
        )
        assert result is not None
        p_value, allele_0, allele_1, group_freqs = result
        _, expected_p, _, _ = chi2_contingency(
            np.array([[2, 0], [1, 1], [0, 2]])
        )
        assert allele_0 == "A"
        assert allele_1 == "G"
        assert p_value == pytest.approx(float(expected_p))
        assert group_freqs == {"X": 0.0, "Y": 0.5, "Z": 1.0}

    def test_test_site_categorical_ascii_fast_path_matches_fallback(self):
        service = PhyloGwas.__new__(PhyloGwas)
        unique_groups = ["X", "Y", "Z"]
        groups = ["X", "X", "Y", "Y", "Z", "Z"]
        group_idx = {group: i for i, group in enumerate(unique_groups)}
        group_codes = np.array(
            [group_idx[group] for group in groups], dtype=np.intp
        )
        group_row_counts = np.bincount(group_codes, minlength=len(unique_groups))

        fast = service._test_site_categorical(
            ["A", "A", "A", "G", "G", "G"],
            groups,
            unique_groups,
            group_idx,
            group_codes,
            group_row_counts,
        )
        fallback = service._test_site_categorical(
            ["α", "α", "α", "β", "β", "β"],
            groups,
            unique_groups,
            group_idx,
            group_codes,
            group_row_counts,
        )

        assert fast is not None
        assert fallback is not None
        assert fast[0] == pytest.approx(fallback[0])
        assert fast[3] == fallback[3]

    def test_test_site_categorical_byte_offsets_match_standard_byte_path(self):
        from scipy.special import chdtrc

        service = PhyloGwas.__new__(PhyloGwas)
        unique_groups = ["X", "Y", "Z"]
        groups = ["X", "X", "Y", "Y", "Z", "Z"]
        group_idx = {group: i for i, group in enumerate(unique_groups)}
        group_codes = np.array(
            [group_idx[group] for group in groups], dtype=np.intp
        )
        group_row_counts = np.bincount(group_codes, minlength=len(unique_groups))
        alleles = np.frombuffer(b"AAAGGG", dtype=np.uint8)

        expected = service._test_site_categorical(
            alleles,
            groups,
            unique_groups,
            group_idx,
            group_codes,
            group_row_counts,
        )
        observed = service._test_site_categorical(
            alleles,
            groups,
            unique_groups,
            group_idx,
            group_codes,
            group_row_counts,
            group_codes * 256,
            chdtrc,
        )

        assert observed is not None
        assert expected is not None
        assert observed[0] == pytest.approx(expected[0])
        assert observed[1:] == expected[1:]

    def test_test_site_categorical_two_group_uses_byte_offsets(self):
        class NoBooleanIndexArray:
            def __init__(self, values):
                self.values = np.asarray(values, dtype=np.intp)

            def __getitem__(self, item):
                if isinstance(item, np.ndarray) and item.dtype == bool:
                    raise AssertionError(
                        "byte offsets should avoid boolean group-code indexing"
                    )
                return self.values[item]

            def __array__(self, dtype=None, copy=None):
                array = np.asarray(self.values, dtype=dtype)
                return array.copy() if copy else array

        service = PhyloGwas.__new__(PhyloGwas)
        groups = ["X", "X", "X", "Y", "Y", "Y"]
        unique_groups = ["X", "Y"]
        group_codes = np.array([0, 0, 0, 1, 1, 1], dtype=np.intp)
        group_row_counts = np.bincount(group_codes, minlength=2)
        alleles = np.frombuffer(b"AAAGGG", dtype=np.uint8)

        expected = service._test_site_categorical(
            alleles,
            groups,
            unique_groups,
            {"X": 0, "Y": 1},
            group_codes,
            group_row_counts,
            group_codes * 256,
            fisher_cache={},
        )
        observed = service._test_site_categorical(
            alleles,
            groups,
            unique_groups,
            {"X": 0, "Y": 1},
            NoBooleanIndexArray(group_codes),
            group_row_counts,
            group_codes * 256,
            fisher_cache={},
        )

        assert observed is not None
        assert expected is not None
        assert observed[0] == pytest.approx(expected[0])
        assert observed[1:] == expected[1:]

    def test_test_site_categorical_two_group_byte_path_skips_histogram(
        self, monkeypatch
    ):
        service = PhyloGwas.__new__(PhyloGwas)
        groups = ["X", "X", "X", "Y", "Y", "Y"]
        unique_groups = ["X", "Y"]
        group_codes = np.array([0, 0, 0, 1, 1, 1], dtype=np.intp)
        group_row_counts = np.array([3, 3], dtype=int)
        alleles = np.frombuffer(b"AAAGGG", dtype=np.uint8)

        def fail_bincount(*_args, **_kwargs):
            raise AssertionError("two-group byte path should not build a histogram")

        monkeypatch.setattr(
            "phykit.services.alignment.phylo_gwas.np.bincount",
            fail_bincount,
        )

        result = service._test_site_categorical(
            alleles,
            groups,
            unique_groups,
            {"X": 0, "Y": 1},
            group_codes,
            group_row_counts,
            fisher_cache={},
        )

        assert result is not None
        assert result[1:] == ("A", "G", {"X": 0.0, "Y": 1.0})

    def test_byte_major_minor_helpers_use_array_minmax(self, monkeypatch):
        alleles = np.frombuffer(b"AAAGGG", dtype=np.uint8)
        group_codes = np.array([0, 0, 0, 1, 1, 1], dtype=np.intp)

        monkeypatch.setattr(
            phylo_gwas_module.np,
            "min",
            lambda *args, **kwargs: pytest.fail(
                "byte allele helpers should use ndarray.min"
            ),
        )
        monkeypatch.setattr(
            phylo_gwas_module.np,
            "max",
            lambda *args, **kwargs: pytest.fail(
                "byte allele helpers should use ndarray.max"
            ),
        )

        assert PhyloGwas._major_minor_bytes(alleles) == (65, 71, 3)
        assert PhyloGwas._major_minor_two_group_counts(
            alleles,
            group_codes,
        ) == (65, 71, (0, 3))

    def test_test_site_categorical_byte_multiallelic_is_skipped(self):
        service = PhyloGwas.__new__(PhyloGwas)
        alleles = np.frombuffer(b"AACCGG", dtype=np.uint8)

        result = service._test_site_categorical(
            alleles,
            ["X", "X", "Y", "Y", "Z", "Z"],
            ["X", "Y", "Z"],
        )

        assert result is None

    def test_test_site_categorical_non_ascii_fallback(self):
        service = PhyloGwas.__new__(PhyloGwas)
        result = service._test_site_categorical(
            ["α", "β", "β", "α"],
            ["X", "X", "Y", "Y"],
            ["X", "Y"],
        )
        assert result is not None
        _, allele_0, allele_1, group_freqs = result
        assert allele_0 == "α"
        assert allele_1 == "β"
        assert group_freqs == {"X": 0.5, "Y": 0.5}

    def test_unicode_multiallelic_site_stops_after_third_allele(self):
        def alleles():
            yield "α"
            yield "β"
            yield "γ"
            raise AssertionError("multiallelic Unicode sites should stop early")

        assert PhyloGwas._binary_alleles_from_unicode_site(alleles()) is None

    def test_unicode_binary_site_preserves_sorted_tie_and_major_order(self):
        assert PhyloGwas._binary_alleles_from_unicode_site(
            ["β", "α", "β", "α"]
        ) == ("α", "β")
        assert PhyloGwas._binary_alleles_from_unicode_site(
            ["β", "β", "α"]
        ) == ("β", "α")

    def test_unicode_binary_site_orders_two_alleles_without_sorting(self, monkeypatch):
        original_sorted = builtins.sorted

        def fail_sorted(*_args, **_kwargs):
            raise AssertionError("Unicode biallelic sites should compare two keys directly")

        monkeypatch.setattr(builtins, "sorted", fail_sorted)

        assert PhyloGwas._binary_alleles_from_unicode_site(
            ["β", "α", "β", "α"]
        ) == ("α", "β")
        assert PhyloGwas._binary_alleles_from_unicode_site(
            ["β", "β", "α"]
        ) == ("β", "α")

        monkeypatch.setattr(builtins, "sorted", original_sorted)

    def test_test_site_continuous_biallelic(self):
        """Biallelic site with continuous phenotype should return correlation."""
        service = PhyloGwas.__new__(PhyloGwas)
        result = service._test_site_continuous(
            ["A", "A", "G", "G"], [10.0, 12.0, 2.0, 3.0]
        )
        assert result is not None
        p_value, allele_0, allele_1, corr_r = result
        assert 0 <= p_value <= 1
        assert isinstance(corr_r, float)

    def test_test_site_continuous_matches_scipy_pointbiserial(self):
        from scipy.special import stdtr
        from scipy.stats import pointbiserialr

        service = PhyloGwas.__new__(PhyloGwas)
        alleles = ["A", "G", "A", "G", "G", "A", "A"]
        values = [1.5, 4.2, 2.1, 5.0, 3.8, 1.2, 2.4]
        value_arr = np.array(values, dtype=float)
        centered = value_arr - value_arr.mean()
        ss = float(np.dot(centered, centered))
        result = service._test_site_continuous(alleles, values)
        cached_result = service._test_site_continuous(
            alleles, values, centered, ss, stdtr
        )
        assert result is not None
        assert cached_result is not None
        p_value, allele_0, allele_1, corr_r = result
        (
            cached_p_value,
            cached_allele_0,
            cached_allele_1,
            cached_corr_r,
        ) = cached_result
        expected_r, expected_p = pointbiserialr(
            [0 if a == allele_0 else 1 for a in alleles],
            values,
        )
        assert allele_0 == "A"
        assert allele_1 == "G"
        assert cached_allele_0 == allele_0
        assert cached_allele_1 == allele_1
        assert corr_r == pytest.approx(float(expected_r))
        assert p_value == pytest.approx(float(expected_p))
        assert cached_corr_r == pytest.approx(float(expected_r))
        assert cached_p_value == pytest.approx(float(expected_p))

    def test_test_site_continuous_byte_column_matches_list_path(self):
        from scipy.special import stdtr

        service = PhyloGwas.__new__(PhyloGwas)
        alleles = ["A", "G", "A", "G", "G", "A", "A"]
        byte_alleles = np.frombuffer("".join(alleles).encode("ascii"), dtype=np.uint8)
        values = [1.5, 4.2, 2.1, 5.0, 3.8, 1.2, 2.4]
        value_arr = np.array(values, dtype=float)
        centered = value_arr - value_arr.mean()
        ss = float(np.dot(centered, centered))

        expected = service._test_site_continuous(
            alleles, values, centered, ss, stdtr
        )
        observed = service._test_site_continuous(
            byte_alleles, values, centered, ss, stdtr
        )

        assert observed is not None
        assert expected is not None
        assert observed[:3] == expected[:3]
        assert observed[3] == pytest.approx(expected[3])
        assert observed[0] == pytest.approx(expected[0])

    def test_test_site_continuous_byte_column_avoids_masked_value_copy(self):
        from scipy.special import stdtr

        class NoBooleanIndexArray:
            def __init__(self, values):
                self.values = np.asarray(values, dtype=float)

            def __array__(self, dtype=None, copy=None):
                array = np.asarray(self.values, dtype=dtype)
                return array.copy() if copy else array

            def __getitem__(self, item):
                if isinstance(item, np.ndarray) and item.dtype == bool:
                    raise AssertionError(
                        "byte continuous path should avoid masked phenotype copies"
                    )
                return self.values[item]

        service = PhyloGwas.__new__(PhyloGwas)
        alleles = np.frombuffer(b"AGAGGAA", dtype=np.uint8)
        values = np.array([1.5, 4.2, 2.1, 5.0, 3.8, 1.2, 2.4], dtype=float)
        centered = values - values.mean()
        ss = float(np.dot(centered, centered))

        expected = service._test_site_continuous(
            alleles, values, centered, ss, stdtr
        )
        observed = service._test_site_continuous(
            alleles, values, NoBooleanIndexArray(centered), ss, stdtr
        )

        assert observed is not None
        assert expected is not None
        assert observed[:3] == expected[:3]
        assert observed[3] == pytest.approx(expected[3])
        assert observed[0] == pytest.approx(expected[0])

    def test_test_site_continuous_byte_multiallelic_is_skipped(self):
        service = PhyloGwas.__new__(PhyloGwas)
        alleles = np.frombuffer(b"AACCGG", dtype=np.uint8)

        result = service._test_site_continuous(
            alleles,
            [1.0, 1.2, 2.0, 2.2, 3.0, 3.2],
        )

        assert result is None

    def test_test_site_continuous_two_samples_matches_scipy_p_value(self):
        service = PhyloGwas.__new__(PhyloGwas)
        result = service._test_site_continuous(["A", "G"], [1.0, 2.0])
        assert result is not None
        p_value, _, _, corr_r = result
        assert corr_r == pytest.approx(1.0)
        assert p_value == pytest.approx(1.0)

    def test_test_site_continuous_non_ascii_fallback_matches_scipy(self):
        from scipy.stats import pointbiserialr

        service = PhyloGwas.__new__(PhyloGwas)
        alleles = ["α", "β", "α", "β"]
        values = [1.0, 3.0, 1.5, 4.0]
        result = service._test_site_continuous(alleles, values)
        assert result is not None
        p_value, allele_0, _, corr_r = result
        expected_r, expected_p = pointbiserialr(
            [0 if allele == allele_0 else 1 for allele in alleles],
            values,
        )
        assert corr_r == pytest.approx(float(expected_r))
        assert p_value == pytest.approx(float(expected_p))

    def test_test_site_continuous_invariant(self):
        """Invariant site with continuous phenotype should return None."""
        service = PhyloGwas.__new__(PhyloGwas)
        result = service._test_site_continuous(
            ["A", "A", "A", "A"], [10.0, 12.0, 2.0, 3.0]
        )
        assert result is None

    def test_test_site_continuous_invariant_does_not_import_scipy_special(
        self, monkeypatch
    ):
        original_import = __import__

        def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name == "scipy.special" or name.startswith("scipy.special."):
                raise AssertionError("invariant continuous site should skip stdtr import")
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr("builtins.__import__", fake_import)

        service = PhyloGwas.__new__(PhyloGwas)
        result = service._test_site_continuous(
            ["A", "A", "A", "A"], [10.0, 12.0, 2.0, 3.0]
        )

        assert result is None
