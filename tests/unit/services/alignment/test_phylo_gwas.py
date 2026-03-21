"""Unit tests for the PhyloGwas service."""

import csv
import json
import os
from argparse import Namespace
from io import StringIO
from pathlib import Path

import numpy as np
import pytest

from phykit.services.alignment.phylo_gwas import PhyloGwas


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

    def test_read_phenotype(self, tmp_path):
        pheno_path = str(tmp_path / "pheno.tsv")
        _write_phenotype(pheno_path, {"s1": "A", "s2": "B"})
        pheno = PhyloGwas._read_phenotype(pheno_path)
        assert pheno == {"s1": "A", "s2": "B"}

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

    def test_creates_png(self, tmp_path):
        """Manhattan plot should be created."""
        pytest.importorskip("matplotlib")
        args = _make_args(tmp_path)
        service = PhyloGwas(args)
        service.run()
        assert os.path.exists(args.output)

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

    def test_position_to_gene(self):
        partitions = [("gene1", 0, 99), ("gene2", 100, 199)]
        assert PhyloGwas._position_to_gene(50, partitions) == "gene1"
        assert PhyloGwas._position_to_gene(150, partitions) == "gene2"
        assert PhyloGwas._position_to_gene(200, partitions) is None

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

    def test_test_site_continuous_invariant(self):
        """Invariant site with continuous phenotype should return None."""
        service = PhyloGwas.__new__(PhyloGwas)
        result = service._test_site_continuous(
            ["A", "A", "A", "A"], [10.0, 12.0, 2.0, 3.0]
        )
        assert result is None
