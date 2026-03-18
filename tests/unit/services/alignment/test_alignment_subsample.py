import os
import pytest
import tempfile
from argparse import Namespace

from phykit.services.alignment.alignment_subsample import AlignmentSubsample
from phykit.errors import PhykitUserError


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _make_gene_list(tmp_path, names=None):
    """Create a gene-list file and return its path."""
    if names is None:
        names = ["gene_a.fa", "gene_b.fa", "gene_c.fa"]
    path = os.path.join(tmp_path, "genes.txt")
    with open(path, "w") as fh:
        for n in names:
            fh.write(f"{n}\n")
    return path


def _make_alignment(tmp_path, sequences=None, fname="aln.fa"):
    """Create a small FASTA alignment and return its path."""
    if sequences is None:
        sequences = {"t1": "ATGATGATG", "t2": "ATGATGATG", "t3": "CTGATGATG"}
    path = os.path.join(tmp_path, fname)
    with open(path, "w") as fh:
        for name, seq in sequences.items():
            fh.write(f">{name}\n{seq}\n")
    return path


def _make_partition(tmp_path, partitions=None, fname="test.partition"):
    """Create a partition file and return its path."""
    if partitions is None:
        partitions = [
            ("gene1", 1, 3),
            ("gene2", 4, 6),
            ("gene3", 7, 9),
        ]
    path = os.path.join(tmp_path, fname)
    with open(path, "w") as fh:
        for name, s, e in partitions:
            fh.write(f"AUTO, {name}={s}-{e}\n")
    return path


def _read_lines(path):
    with open(path) as fh:
        return [l.strip() for l in fh if l.strip()]


def _read_fasta(path):
    seqs = {}
    current = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                current = line[1:]
                seqs[current] = ""
            elif current is not None:
                seqs[current] += line
    return seqs


# ---------------------------------------------------------------------------
# genes mode
# ---------------------------------------------------------------------------

class TestGenesMode:
    def test_genes_mode_number(self, tmp_path):
        gene_list = _make_gene_list(str(tmp_path))
        prefix = os.path.join(str(tmp_path), "out")
        args = Namespace(
            mode="genes", alignment=None, list=gene_list, partition=None,
            number=2, fraction=None, seed=42, bootstrap=False,
            output=prefix, json=False,
        )
        AlignmentSubsample(args).run()
        lines = _read_lines(f"{prefix}.txt")
        assert len(lines) == 2

    def test_genes_mode_fraction(self, tmp_path):
        gene_list = _make_gene_list(str(tmp_path))
        prefix = os.path.join(str(tmp_path), "out")
        args = Namespace(
            mode="genes", alignment=None, list=gene_list, partition=None,
            number=None, fraction=0.67, seed=42, bootstrap=False,
            output=prefix, json=False,
        )
        AlignmentSubsample(args).run()
        lines = _read_lines(f"{prefix}.txt")
        # round(0.67 * 3) = 2
        assert len(lines) == 2

    def test_genes_mode_seed_reproducible(self, tmp_path):
        gene_list = _make_gene_list(str(tmp_path))
        prefix1 = os.path.join(str(tmp_path), "out1")
        prefix2 = os.path.join(str(tmp_path), "out2")
        for prefix in (prefix1, prefix2):
            args = Namespace(
                mode="genes", alignment=None, list=gene_list, partition=None,
                number=2, fraction=None, seed=123, bootstrap=False,
                output=prefix, json=False,
            )
            AlignmentSubsample(args).run()
        assert _read_lines(f"{prefix1}.txt") == _read_lines(f"{prefix2}.txt")

    def test_genes_mode_bootstrap(self, tmp_path):
        """With replacement, selecting N >= total can include duplicates."""
        gene_list = _make_gene_list(str(tmp_path))
        prefix = os.path.join(str(tmp_path), "out")
        args = Namespace(
            mode="genes", alignment=None, list=gene_list, partition=None,
            number=6, fraction=None, seed=42, bootstrap=True,
            output=prefix, json=False,
        )
        AlignmentSubsample(args).run()
        lines = _read_lines(f"{prefix}.txt")
        assert len(lines) == 6
        # Some should be duplicates since we pick 6 from 3
        assert len(set(lines)) < len(lines)


# ---------------------------------------------------------------------------
# sites mode
# ---------------------------------------------------------------------------

class TestSitesMode:
    def test_sites_mode_number(self, tmp_path):
        aln = _make_alignment(str(tmp_path))
        prefix = os.path.join(str(tmp_path), "out")
        args = Namespace(
            mode="sites", alignment=aln, list=None, partition=None,
            number=5, fraction=None, seed=42, bootstrap=False,
            output=prefix, json=False,
        )
        AlignmentSubsample(args).run()
        seqs = _read_fasta(f"{prefix}.fa")
        for seq in seqs.values():
            assert len(seq) == 5

    def test_sites_mode_fraction(self, tmp_path):
        aln = _make_alignment(str(tmp_path))
        prefix = os.path.join(str(tmp_path), "out")
        args = Namespace(
            mode="sites", alignment=aln, list=None, partition=None,
            number=None, fraction=0.5, seed=42, bootstrap=False,
            output=prefix, json=False,
        )
        AlignmentSubsample(args).run()
        seqs = _read_fasta(f"{prefix}.fa")
        # round(0.5 * 9) = 4 or 5
        expected = round(0.5 * 9)
        for seq in seqs.values():
            assert len(seq) == expected

    def test_sites_mode_seed_reproducible(self, tmp_path):
        aln = _make_alignment(str(tmp_path))
        prefix1 = os.path.join(str(tmp_path), "out1")
        prefix2 = os.path.join(str(tmp_path), "out2")
        for prefix in (prefix1, prefix2):
            args = Namespace(
                mode="sites", alignment=aln, list=None, partition=None,
                number=5, fraction=None, seed=99, bootstrap=False,
                output=prefix, json=False,
            )
            AlignmentSubsample(args).run()
        assert _read_fasta(f"{prefix1}.fa") == _read_fasta(f"{prefix2}.fa")


# ---------------------------------------------------------------------------
# partitions mode
# ---------------------------------------------------------------------------

class TestPartitionsMode:
    def test_partitions_mode_number(self, tmp_path):
        seqs = {"t1": "ATGATGATG", "t2": "CTGATGCTG", "t3": "ATGATCATG"}
        aln = _make_alignment(str(tmp_path), seqs)
        part = _make_partition(str(tmp_path))
        prefix = os.path.join(str(tmp_path), "out")
        args = Namespace(
            mode="partitions", alignment=aln, list=None, partition=part,
            number=2, fraction=None, seed=42, bootstrap=False,
            output=prefix, json=False,
        )
        AlignmentSubsample(args).run()
        seqs_out = _read_fasta(f"{prefix}.fa")
        # 2 partitions of width 3 each => length 6
        for seq in seqs_out.values():
            assert len(seq) == 6

    def test_partitions_mode_renumbers(self, tmp_path):
        seqs = {"t1": "ATGATGATG", "t2": "CTGATGCTG"}
        aln = _make_alignment(str(tmp_path), seqs)
        part = _make_partition(str(tmp_path))
        prefix = os.path.join(str(tmp_path), "out")
        args = Namespace(
            mode="partitions", alignment=aln, list=None, partition=part,
            number=2, fraction=None, seed=42, bootstrap=False,
            output=prefix, json=False,
        )
        AlignmentSubsample(args).run()
        part_lines = _read_lines(f"{prefix}.partition")
        assert len(part_lines) == 2
        # First partition should start at 1
        assert "=1-3" in part_lines[0]
        # Second partition should start at 4
        assert "=4-6" in part_lines[1]

    def test_partitions_bootstrap_duplicates(self, tmp_path):
        seqs = {"t1": "ATGATGATG", "t2": "CTGATGCTG"}
        aln = _make_alignment(str(tmp_path), seqs)
        part = _make_partition(str(tmp_path))
        prefix = os.path.join(str(tmp_path), "out")
        args = Namespace(
            mode="partitions", alignment=aln, list=None, partition=part,
            number=6, fraction=None, seed=42, bootstrap=True,
            output=prefix, json=False,
        )
        AlignmentSubsample(args).run()
        seqs_out = _read_fasta(f"{prefix}.fa")
        # 6 partitions of 3 columns = 18
        for seq in seqs_out.values():
            assert len(seq) == 18
        part_lines = _read_lines(f"{prefix}.partition")
        assert len(part_lines) == 6


# ---------------------------------------------------------------------------
# validation errors
# ---------------------------------------------------------------------------

class TestValidation:
    def test_number_and_fraction_exclusive(self, tmp_path):
        with pytest.raises(PhykitUserError):
            args = Namespace(
                mode="genes", alignment=None, list="x.txt", partition=None,
                number=2, fraction=0.5, seed=None, bootstrap=False,
                output="out", json=False,
            )
            AlignmentSubsample(args)

    def test_neither_number_nor_fraction(self, tmp_path):
        with pytest.raises(PhykitUserError):
            args = Namespace(
                mode="genes", alignment=None, list="x.txt", partition=None,
                number=None, fraction=None, seed=None, bootstrap=False,
                output="out", json=False,
            )
            AlignmentSubsample(args)

    def test_genes_requires_list(self, tmp_path):
        prefix = os.path.join(str(tmp_path), "out")
        args = Namespace(
            mode="genes", alignment=None, list=None, partition=None,
            number=1, fraction=None, seed=42, bootstrap=False,
            output=prefix, json=False,
        )
        with pytest.raises(PhykitUserError):
            AlignmentSubsample(args).run()

    def test_partitions_requires_alignment_and_partition(self, tmp_path):
        prefix = os.path.join(str(tmp_path), "out")
        args = Namespace(
            mode="partitions", alignment=None, list=None, partition=None,
            number=1, fraction=None, seed=42, bootstrap=False,
            output=prefix, json=False,
        )
        with pytest.raises(PhykitUserError):
            AlignmentSubsample(args).run()

    def test_sites_requires_alignment(self, tmp_path):
        prefix = os.path.join(str(tmp_path), "out")
        args = Namespace(
            mode="sites", alignment=None, list=None, partition=None,
            number=1, fraction=None, seed=42, bootstrap=False,
            output=prefix, json=False,
        )
        with pytest.raises(PhykitUserError):
            AlignmentSubsample(args).run()
