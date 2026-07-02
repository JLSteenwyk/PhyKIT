import os
import pytest
import subprocess
import sys
import tempfile
from argparse import Namespace
from operator import itemgetter

from phykit.services.alignment.alignment_subsample import AlignmentSubsample
from phykit.errors import PhykitUserError


def test_module_import_does_not_import_biopython_fasta_parser():
    code = """
import sys
import phykit.services.alignment.alignment_subsample
assert "typing" not in sys.modules
assert "random" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "Bio.SeqIO.FastaIO" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


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


class TestPartitionParsing:
    def test_parse_partition_file_handles_comments_whitespace_and_trailing_text(self, tmp_path):
        path = os.path.join(str(tmp_path), "test.partition")
        with open(path, "w") as fh:
            fh.write("   # ignored\n")
            fh.write("\n")
            fh.write("AUTO, gene1=1-3\n")
            fh.write("DNA,   gene2   =   4   -   6 trailing text\n")
            fh.write("invalid line\n")
            fh.write("DNA model, gene3=7-9\n")
            fh.write("AUTO, gene 4=10-12\n")

        assert AlignmentSubsample._parse_partition_file(path) == [
            ("gene1", 1, 3),
            ("gene2", 4, 6),
        ]

    def test_read_list_file_skips_whitespace_prefixed_comments(self, tmp_path):
        path = os.path.join(str(tmp_path), "taxa.txt")
        with open(path, "w") as fh:
            fh.write("   # ignored\n")
            fh.write("\n")
            fh.write("taxon_a\n")
            fh.write("\t# ignored\n")
            fh.write("taxon_b\n")

        assert AlignmentSubsample._read_list_file(path) == ["taxon_a", "taxon_b"]


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

    def test_genes_mode_writes_selected_paths_once_with_trailing_newline(
        self, tmp_path, monkeypatch
    ):
        gene_list = _make_gene_list(str(tmp_path))
        prefix = os.path.join(str(tmp_path), "out")
        args = Namespace(
            mode="genes", alignment=None, list=gene_list, partition=None,
            number=2, fraction=None, seed=123, bootstrap=False,
            output=prefix, json=False,
        )
        service = AlignmentSubsample(args)
        monkeypatch.setattr(
            service,
            "_sample",
            lambda _rng, _items, _n: ["gene3.fa", "gene1.fa"],
        )

        service.run()

        assert open(f"{prefix}.txt").read() == "gene3.fa\ngene1.fa\n"

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
    def test_select_sites_single_site(self):
        assert AlignmentSubsample._select_sites("ACGT", itemgetter(2)) == "G"

    def test_select_sites_preserves_repeated_indices(self):
        assert AlignmentSubsample._select_sites("ACGT", itemgetter(2, 2, 0)) == "GGA"

    def test_selected_index_ranges_groups_contiguous_sites(self):
        assert AlignmentSubsample._selected_index_ranges(
            [0, 1, 2, 5, 7, 8]
        ) == [(0, 3), (5, 6), (7, 9)]

    def test_select_site_ranges_joins_sequence_slices(self):
        assert AlignmentSubsample._select_site_ranges(
            "ACGTACGT",
            [(0, 2), (4, 8)],
        ) == "ACACGT"

    def test_select_site_ranges_handles_many_short_ranges(self):
        sequence = "ACGT" * 100
        ranges = [(idx, idx + 1) for idx in range(0, len(sequence), 2)]
        assert (
            AlignmentSubsample._select_site_ranges(sequence, ranges)
            == sequence[::2]
        )

    def test_read_alignment_uses_first_header_token_and_last_duplicate(self, tmp_path):
        aln = os.path.join(str(tmp_path), "aln.fa")
        with open(aln, "w") as fh:
            fh.write(">t1 first description\nAA AA\n")
            fh.write(">t2 second description\nCc c\nC\n")
            fh.write(">t1 replacement description\nGGGG\n")

        assert AlignmentSubsample._read_alignment(aln) == {
            "t1": "GGGG",
            "t2": "CccC",
        }

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

    def test_sites_mode_full_nonbootstrap_reuses_parsed_sequences(
        self, tmp_path, monkeypatch
    ):
        aln = _make_alignment(str(tmp_path))
        prefix = os.path.join(str(tmp_path), "out")
        args = Namespace(
            mode="sites", alignment=aln, list=None, partition=None,
            number=9, fraction=None, seed=42, bootstrap=False,
            output=prefix, json=False,
        )
        monkeypatch.setattr(
            AlignmentSubsample,
            "_select_sites",
            staticmethod(
                lambda *_args, **_kwargs: (_ for _ in ()).throw(
                    AssertionError("full-site selection should not rebuild sequences")
                )
            ),
        )

        AlignmentSubsample(args).run()

        assert _read_fasta(f"{prefix}.fa") == _read_fasta(aln)

    def test_sites_mode_full_nonbootstrap_writes_original_sequence_mapping(
        self, tmp_path, monkeypatch
    ):
        sequences = {"t1": "AAAA", "t2": "CCCC"}
        captured = {}
        args = Namespace(
            mode="sites", alignment="aln.fa", list=None, partition=None,
            number=4, fraction=None, seed=42, bootstrap=False,
            output=os.path.join(str(tmp_path), "out"), json=False,
        )
        service = AlignmentSubsample(args)
        monkeypatch.setattr(service, "_read_alignment", lambda _path: sequences)
        monkeypatch.setattr(
            service,
            "_write_fasta",
            lambda _path, written_sequences: captured.setdefault(
                "sequences", written_sequences
            ),
        )
        monkeypatch.setattr(service, "_print_summary", lambda *_args: None)

        service._run_sites(None)

        assert captured["sequences"] is sequences

    def test_sites_mode_clustered_nonbootstrap_selection_uses_ranges(
        self, tmp_path, monkeypatch
    ):
        class FakeRng:
            def sample(self, _items, k):
                assert k == 9
                return [8, 7, 6, 5, 4, 3, 2, 1, 0]

        captured = {}
        args = Namespace(
            mode="sites", alignment="aln.fa", list=None, partition=None,
            number=9, fraction=None, seed=42, bootstrap=False,
            output=os.path.join(str(tmp_path), "out"), json=False,
        )
        service = AlignmentSubsample(args)
        monkeypatch.setattr(
            service,
            "_read_alignment",
            lambda _path: {"t1": "ACGTACGTAA", "t2": "TTTTCCCCGG"},
        )
        monkeypatch.setattr(
            AlignmentSubsample,
            "_select_sites",
            staticmethod(
                lambda *_args, **_kwargs: (_ for _ in ()).throw(
                    AssertionError("clustered site selection should use ranges")
                )
            ),
        )
        monkeypatch.setattr(
            service,
            "_write_fasta",
            lambda _path, written_sequences: captured.setdefault(
                "sequences", written_sequences
            ),
        )
        monkeypatch.setattr(service, "_print_summary", lambda *_args: None)

        service._run_sites(FakeRng())

        assert captured["sequences"] == {"t1": "ACGTACGTA", "t2": "TTTTCCCCG"}

    def test_sites_mode_stops_at_first_length_mismatch(
        self, tmp_path, monkeypatch
    ):
        class EarlyMismatchSequences(dict):
            def values(self):
                yield "AAAA"
                yield "AAA"
                raise AssertionError(
                    "length validation should stop at the first mismatch"
                )

        args = Namespace(
            mode="sites", alignment="aln.fa", list=None, partition=None,
            number=2, fraction=None, seed=42, bootstrap=False,
            output=os.path.join(str(tmp_path), "out"), json=False,
        )
        service = AlignmentSubsample(args)
        monkeypatch.setattr(
            service,
            "_read_alignment",
            lambda _path: EarlyMismatchSequences(
                {"t1": "AAAA", "t2": "AAA", "t3": "AAAA"}
            ),
        )

        with pytest.raises(PhykitUserError) as exc_info:
            service._run_sites(None)

        assert exc_info.value.messages == [
            "All sequences in the alignment must have the same length."
        ]

    def test_write_fasta_batches_records_preserving_exact_text(
        self, tmp_path, monkeypatch
    ):
        monkeypatch.setattr(
            "phykit.services.alignment.alignment_subsample._FASTA_WRITE_CHUNK_ROWS",
            2,
        )
        path = os.path.join(str(tmp_path), "out.fa")

        AlignmentSubsample._write_fasta(
            path,
            {"t1": "ACGT", "t2": "TGCA", "t3": "NNNN"},
        )

        assert open(path).read() == ">t1\nACGT\n>t2\nTGCA\n>t3\nNNNN\n"


# ---------------------------------------------------------------------------
# partitions mode
# ---------------------------------------------------------------------------

class TestPartitionsMode:
    def test_assemble_partition_subsample_preserves_selection_and_duplicates(self):
        sequences = {"t1": "AAACCCGGG", "t2": "TTTGGGCCC"}
        selected = [
            ("gene3", 7, 9),
            ("gene1", 1, 3),
            ("gene3", 7, 9),
        ]

        new_sequences, new_partitions = AlignmentSubsample._assemble_partition_subsample(
            sequences,
            selected,
        )

        assert new_sequences == {"t1": "GGGAAAGGG", "t2": "CCCTTTCCC"}
        assert new_partitions == [
            ("gene3", 1, 3),
            ("gene1", 4, 6),
            ("gene3_dup1", 7, 9),
        ]

    def test_write_partition_file_batches_rows_preserving_exact_text(
        self, tmp_path, monkeypatch
    ):
        monkeypatch.setattr(
            "phykit.services.alignment.alignment_subsample._PARTITION_WRITE_CHUNK_ROWS",
            2,
        )
        path = os.path.join(str(tmp_path), "out.partition")

        AlignmentSubsample._write_partition_file(
            path,
            [("gene1", 1, 3), ("gene2", 4, 6), ("gene3", 7, 9)],
        )

        assert (
            open(path).read()
            == "AUTO, gene1=1-3\nAUTO, gene2=4-6\nAUTO, gene3=7-9\n"
        )

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

    def test_partitions_bootstrap_repeats_selected_sequence_slices(self, tmp_path):
        seqs = {"t1": "AAACCCGGG", "t2": "TTTGGGCCC"}
        aln = _make_alignment(str(tmp_path), seqs)
        part = _make_partition(str(tmp_path))
        prefix = os.path.join(str(tmp_path), "out")
        args = Namespace(
            mode="partitions", alignment=aln, list=None, partition=part,
            number=3, fraction=None, seed=1, bootstrap=True,
            output=prefix, json=False,
        )
        AlignmentSubsample(args).run()
        seqs_out = _read_fasta(f"{prefix}.fa")
        part_lines = _read_lines(f"{prefix}.partition")

        assert seqs_out["t1"] == "AAAGGGGGG"
        assert seqs_out["t2"] == "TTTCCCCCC"
        assert part_lines == [
            "AUTO, gene1=1-3",
            "AUTO, gene3=4-6",
            "AUTO, gene3_dup1=7-9",
        ]


class TestSummaryOutput:
    def test_print_summary_batches_text_output(self, monkeypatch):
        svc = AlignmentSubsample.__new__(AlignmentSubsample)
        svc.seed = None
        svc.bootstrap = True
        svc.json_output = False
        printed = []

        def fake_print(*args, **kwargs):
            printed.append((args, kwargs))

        monkeypatch.setattr("builtins.print", fake_print)

        svc._print_summary(
            "partitions",
            5,
            3,
            ["out.fa", "out.partition"],
        )

        expected = "\n".join([
            "Alignment Subsampling",
            "Mode: partitions",
            "Total partitions: 5",
            "Selected: 3",
            "Bootstrap: yes",
            "Seed: random",
            "Output: out.fa",
            "Output: out.partition",
        ])
        assert printed == [((expected,), {})]


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
