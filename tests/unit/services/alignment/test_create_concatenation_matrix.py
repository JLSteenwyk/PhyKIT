from argparse import Namespace
from collections import defaultdict
from pathlib import Path
import subprocess
import sys

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import pytest

from phykit.services.alignment.create_concatenation_matrix import (
    CreateConcatenationMatrix,
    _ParsedFastaRecord,
)
import phykit.services.alignment.create_concatenation_matrix as ccm_module


def test_module_import_does_not_import_biopython_fasta_parser():
    code = """
import sys
import phykit.services.alignment.create_concatenation_matrix as module
assert hasattr(module.np, "__getattr__")
assert hasattr(module.mp, "cpu_count")
assert callable(module.read_single_column_file_to_list)
assert module._OCCUPANCY_STATE_LOOKUP is None
assert "typing" not in sys.modules
assert "numpy" not in sys.modules
assert "concurrent.futures" not in sys.modules
assert "multiprocessing" not in sys.modules
assert "json" not in sys.modules
assert "textwrap" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
assert "phykit.helpers.files" not in sys.modules
assert "Bio.SeqIO.FastaIO" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def _write_fasta(path: Path, records):
    lines = []
    for taxon, seq in records:
        lines.append(f">{taxon}")
        lines.append(seq)
    path.write_text("\n".join(lines) + "\n")


@pytest.fixture
def args():
    return Namespace(alignment_list="/some/path/to/file", prefix="some_prefix")


class TestCreateConcatenationMatrix:
    def test_init_sets_alignment_list_path(self, args):
        ccm = CreateConcatenationMatrix(args)
        assert ccm.alignment_list_path == args.alignment_list
        assert ccm.prefix == args.prefix
        assert ccm.output_file_path is None
        assert ccm.json_output is False
        assert ccm.plot_occupancy is False
        assert ccm.plot_output is None
        assert ccm.threshold == 0
        assert ccm.plot_config is not None

    def test_process_args_reads_optional_flags(self):
        parsed = CreateConcatenationMatrix(
            Namespace(
                alignment_list="x.list",
                prefix="x",
                json=True,
                plot_occupancy=True,
                plot_output="occ.png",
            )
        ).process_args(
            Namespace(
                alignment_list="x.list",
                prefix="x",
                json=True,
                plot_occupancy=True,
                plot_output="occ.png",
            )
        )
        assert parsed["json_output"] is True
        assert parsed["plot_occupancy"] is True
        assert parsed["plot_output"] == "occ.png"

    def test_read_alignment_paths_missing_file_exits(self, args):
        ccm = CreateConcatenationMatrix(args)
        with pytest.raises(SystemExit) as exc:
            ccm.read_alignment_paths("/definitely/not/found.txt")
        assert exc.value.code == 2

    def test_get_taxa_from_alignment(self, tmp_path, args):
        fasta = tmp_path / "gene.fa"
        _write_fasta(fasta, [("A", "ACGT"), ("B", "A-GT")])
        ccm = CreateConcatenationMatrix(args)
        taxa = ccm._get_taxa_from_alignment(str(fasta))
        assert taxa == {"A", "B"}

    def test_get_taxa_from_alignment_uses_first_header_token(self, tmp_path, args):
        fasta = tmp_path / "gene.fa"
        fasta.write_text(
            ">A description text\nAC GT\n>B another description\nA- GT\n"
        )
        ccm = CreateConcatenationMatrix(args)
        taxa = ccm._get_taxa_from_alignment(str(fasta))
        assert taxa == {"A", "B"}

    def test_get_taxa_names_sequential(self, tmp_path, args):
        fastas = []
        for idx, entries in enumerate(
            [
                [("A", "AAAA"), ("B", "CCCC")],
                [("A", "GGGG"), ("C", "TTTT")],
            ]
        ):
            fasta = tmp_path / f"gene_{idx}.fa"
            _write_fasta(fasta, entries)
            fastas.append(str(fasta))

        ccm = CreateConcatenationMatrix(args)
        taxa = ccm.get_taxa_names(fastas)
        assert taxa == ["A", "B", "C"]

    def test_get_taxa_names_parallel_fallback(self, monkeypatch, args):
        class FailingExecutor:
            def __init__(self, *_, **__):
                pass

            def __enter__(self):
                raise RuntimeError("no multiprocessing")

            def __exit__(self, exc_type, exc, tb):
                return False

        ccm = CreateConcatenationMatrix(args)
        monkeypatch.setattr(ccm_module, "ProcessPoolExecutor", FailingExecutor)
        monkeypatch.setattr(
            ccm,
            "_get_taxa_from_alignment",
            lambda path: {path.split("_")[-1]},
        )

        paths = [f"gene_{i}" for i in range(11)]
        taxa = ccm.get_taxa_names(paths)
        assert len(taxa) == 11

    def test_get_list_of_taxa_and_records_uses_first_header_token_and_keeps_duplicates(
        self, tmp_path, args
    ):
        fasta = tmp_path / "gene.fa"
        fasta.write_text(
            ">A description\nACGT\n"
            ">A duplicate\nTT TT\nAA\n"
            ">B other text\nCC\nCC\n"
        )
        ccm = CreateConcatenationMatrix(args)

        taxa, records = ccm.get_list_of_taxa_and_records(str(fasta))

        assert taxa == {"A", "B"}
        assert [record.id for record in records] == ["A", "A", "B"]
        assert [str(record.seq) for record in records] == ["ACGT", "TTTTAA", "CCCC"]
        assert not hasattr(records[0], "__dict__")

    def test_create_missing_seq_str(self, args):
        ccm = CreateConcatenationMatrix(args)
        missing_seq, og_len = ccm.create_missing_seq_str([SeqRecord(Seq("ACGT"), id="A")])
        assert missing_seq == "????"
        assert og_len == 4

    def test_create_missing_seq_str_empty_exits(self, args):
        ccm = CreateConcatenationMatrix(args)
        with pytest.raises(SystemExit) as exc:
            ccm.create_missing_seq_str([])
        assert exc.value.code == 2

    def test_process_taxa_sequences_adds_missing(self, args):
        ccm = CreateConcatenationMatrix(args)
        records = [SeqRecord(Seq("ACGT"), id="A")]
        concatenated_seqs = defaultdict(list)
        ccm.process_taxa_sequences(records, ["A", "B"], concatenated_seqs, "????")
        assert concatenated_seqs["A"] == ["ACGT"]
        assert concatenated_seqs["B"] == ["????"]

    def test_process_taxa_sequences_uses_cached_taxa_set(self, args):
        ccm = CreateConcatenationMatrix(args)
        records = [SeqRecord(Seq("ACGT"), id="A")]
        concatenated_seqs = defaultdict(list)

        ccm.process_taxa_sequences(
            records,
            ["A", "B"],
            concatenated_seqs,
            "????",
            taxa_set={"A", "B"},
        )

        assert concatenated_seqs["A"] == ["ACGT"]
        assert concatenated_seqs["B"] == ["????"]

    def test_process_taxa_sequences_reuses_present_taxa(self, args):
        class NoAddSet(set):
            def add(self, item):
                raise AssertionError("present taxa should not be rebuilt")

        ccm = CreateConcatenationMatrix(args)
        records = [SeqRecord(Seq("ACGT"), id="A")]
        concatenated_seqs = defaultdict(list)

        ccm.process_taxa_sequences(
            records,
            ["A", "B"],
            concatenated_seqs,
            "????",
            taxa_set={"A", "B"},
            present_taxa=NoAddSet({"A"}),
        )

        assert concatenated_seqs["A"] == ["ACGT"]
        assert concatenated_seqs["B"] == ["????"]

    def test_process_taxa_sequences_all_present_skips_missing_scan(self, args):
        class NoIterTaxa(list):
            def __iter__(self):
                raise AssertionError("complete present taxa should skip missing scan")

        ccm = CreateConcatenationMatrix(args)
        records = [
            SeqRecord(Seq("ACGT"), id="A"),
            SeqRecord(Seq("TGCA"), id="B"),
        ]
        concatenated_seqs = defaultdict(list)

        ccm.process_taxa_sequences(
            records,
            NoIterTaxa(["A", "B"]),
            concatenated_seqs,
            "????",
            taxa_set={"A", "B"},
            present_taxa={"A", "B"},
        )

        assert concatenated_seqs["A"] == ["ACGT"]
        assert concatenated_seqs["B"] == ["TGCA"]

    def test_process_taxa_sequences_preserves_duplicate_records(self, args):
        ccm = CreateConcatenationMatrix(args)
        records = [
            SeqRecord(Seq("ACGT"), id="A"),
            SeqRecord(Seq("TGCA"), id="A"),
        ]
        concatenated_seqs = defaultdict(list)

        ccm.process_taxa_sequences(
            records,
            ["A", "B"],
            concatenated_seqs,
            "????",
            taxa_set={"A", "B"},
        )

        assert concatenated_seqs["A"] == ["ACGT", "TGCA"]
        assert concatenated_seqs["B"] == ["????"]

    def test_process_taxa_sequences_converts_mixed_external_records(self, args):
        class Record:
            def __init__(self, record_id, sequence):
                self.id = record_id
                self.seq = sequence

        ccm = CreateConcatenationMatrix(args)
        records = [
            Record("A", "ACGT"),
            Record("B", Seq("TGCA")),
        ]
        concatenated_seqs = defaultdict(list)

        ccm.process_taxa_sequences(
            records,
            ["A", "B"],
            concatenated_seqs,
            "????",
            taxa_set={"A", "B"},
        )

        assert concatenated_seqs["A"] == ["ACGT"]
        assert concatenated_seqs["B"] == ["TGCA"]

    def test_process_taxa_sequences_accepts_parsed_fasta_records(self, args):
        ccm = CreateConcatenationMatrix(args)
        records = [_ParsedFastaRecord("A", "ACGT")]
        concatenated_seqs = defaultdict(list)

        ccm.process_taxa_sequences(
            records,
            ["A", "B"],
            concatenated_seqs,
            "????",
            taxa_set={"A", "B"},
        )

        assert concatenated_seqs["A"] == ["ACGT"]
        assert concatenated_seqs["B"] == ["????"]

    def test_process_taxa_sequences_adds_missing_in_taxa_order(self, args):
        ccm = CreateConcatenationMatrix(args)
        records = [SeqRecord(Seq("ACGT"), id="B")]
        concatenated_seqs = defaultdict(list)

        ccm.process_taxa_sequences(
            records,
            ["A", "B", "C"],
            concatenated_seqs,
            "????",
            taxa_set={"A", "B", "C"},
        )

        assert list(concatenated_seqs.keys()) == ["B", "A", "C"]
        assert concatenated_seqs["A"] == ["????"]
        assert concatenated_seqs["B"] == ["ACGT"]
        assert concatenated_seqs["C"] == ["????"]

    def test_process_taxa_sequences_duplicate_taxa_list_adds_missing_once(self, args):
        ccm = CreateConcatenationMatrix(args)
        records = [SeqRecord(Seq("ACGT"), id="A")]
        concatenated_seqs = defaultdict(list)

        ccm.process_taxa_sequences(
            records,
            ["A", "B", "B"],
            concatenated_seqs,
            "????",
            taxa_set={"A", "B"},
        )

        assert concatenated_seqs["A"] == ["ACGT"]
        assert concatenated_seqs["B"] == ["????"]

    def test_process_taxa_sequences_uses_precomputed_missing_taxa(self, args):
        class NoIterTaxa(list):
            def __iter__(self):
                raise AssertionError("precomputed missing taxa should avoid taxa scan")

        ccm = CreateConcatenationMatrix(args)
        records = [_ParsedFastaRecord("A", "ACGT")]
        concatenated_seqs = defaultdict(list)

        ccm.process_taxa_sequences(
            records,
            NoIterTaxa(["A", "B", "C"]),
            concatenated_seqs,
            "????",
            taxa_set={"A", "B", "C"},
            present_taxa={"A"},
            missing_taxa=["B", "C"],
        )

        assert concatenated_seqs["A"] == ["ACGT"]
        assert concatenated_seqs["B"] == ["????"]
        assert concatenated_seqs["C"] == ["????"]

    def test_add_to_partition_info(self, args):
        ccm = CreateConcatenationMatrix(args)
        partition_info, first_len, second_len = ccm.add_to_partition_info(
            [], 10, "AUTO", "g1.fa", 1, 0
        )
        assert partition_info == ["AUTO, g1.fa=1-10\n"]
        assert first_len == 11
        assert second_len == 10

    def test_add_to_occupancy_info(self, args):
        ccm = CreateConcatenationMatrix(args)
        occ = ccm.add_to_occupancy_info([], {"A"}, ["A", "B"], "g1.fa")
        assert occ == ["g1.fa\t1\t1\t0.5000\tB\n"]

    def test_add_to_occupancy_info_preserves_sorted_unique_missing_taxa(self, args):
        ccm = CreateConcatenationMatrix(args)

        occ = ccm.add_to_occupancy_info([], {"A"}, ["B", "A", "B", "C"], "g1.fa")

        assert occ == ["g1.fa\t1\t2\t0.2500\tB;C\n"]

    def test_add_to_occupancy_info_uses_cached_sorted_taxa(self, args):
        ccm = CreateConcatenationMatrix(args)

        occ = ccm.add_to_occupancy_info(
            [],
            {"A", "C"},
            ["A", "B", "C", "D"],
            "g1.fa",
            sorted_taxa=["A", "B", "C", "D"],
            total_taxa_count=4,
        )

        assert occ == ["g1.fa\t2\t2\t0.5000\tB;D\n"]

    def test_add_to_occupancy_info_uses_precomputed_missing_taxa(self, args):
        class NoIterTaxa(list):
            def __iter__(self):
                raise AssertionError("precomputed missing taxa should avoid taxa scan")

        ccm = CreateConcatenationMatrix(args)

        occ = ccm.add_to_occupancy_info(
            [],
            {"A", "C"},
            NoIterTaxa(["A", "B", "C", "D"]),
            "g1.fa",
            total_taxa_count=4,
            missing_taxa=["B", "D"],
        )

        assert occ == ["g1.fa\t2\t2\t0.5000\tB;D\n"]

    def test_fasta_and_text_file_writers(self, tmp_path, args):
        ccm = CreateConcatenationMatrix(args)
        fasta_out = tmp_path / "out.fa"
        ccm.fasta_file_write(str(fasta_out), {"A": ["AA", "CC"], "B": ["GG"]})
        assert fasta_out.read_text() == ">A\nAACC\n>B\nGG\n"

        text_out = tmp_path / "out.txt"
        ccm.write_occupancy_or_partition_file(["x\n", "y\n"], str(text_out))
        assert text_out.read_text() == "x\ny\n"

    def test_fasta_file_write_uses_chunked_writes(self, monkeypatch, args):
        ccm = CreateConcatenationMatrix(args)

        class TrackingHandle:
            def __init__(self):
                self.chunks = []

            def __enter__(self):
                return self

            def __exit__(self, *_args):
                return False

            def write(self, text):
                self.chunks.append(text)

            def writelines(self, rows):
                raise AssertionError("FASTA output should use chunked writes")

        handle = TrackingHandle()

        def fake_open(path, mode, buffering=None):
            assert path == "out.fa"
            assert mode == "w"
            assert buffering == 8192
            return handle

        monkeypatch.setattr(ccm_module, "open", fake_open, raising=False)
        monkeypatch.setattr(ccm_module, "_FASTA_WRITE_CHUNK_ROWS", 1)

        ccm.fasta_file_write("out.fa", {"A": ["AA", "CC"], "B": ["GG"]})

        assert handle.chunks == [">A\nAACC\n", ">B\nGG\n"]

    def test_append_ordered_sequences_uses_cached_taxon_lists(self, args):
        ccm = CreateConcatenationMatrix(args)
        taxon_seq_lists = [[], []]

        ccm._append_ordered_sequences(
            ["B", "A"],
            {"A": "AAAA", "B": "BBBB"},
            taxon_seq_lists,
        )

        assert taxon_seq_lists == [["BBBB"], ["AAAA"]]

    def test_process_alignment_file_populates_missing_taxa(self, tmp_path):
        fasta = tmp_path / "gene.fa"
        _write_fasta(fasta, [("A", "ACGT"), ("C", "AAGT")])
        _, seq_dict, present_taxa, og_len = CreateConcatenationMatrix._process_alignment_file(
            str(fasta),
            ["A", "B", "C"],
        )
        assert present_taxa == {"A", "C"}
        assert og_len == 4
        assert seq_dict["A"] == "ACGT"
        assert seq_dict["B"] == "????"
        assert seq_dict["C"] == "AAGT"

    def test_process_alignment_file_preserves_first_duplicate_record(self, tmp_path):
        fasta = tmp_path / "gene.fa"
        _write_fasta(fasta, [("A", "ACGT"), ("A", "TTTT"), ("B", "CCCC")])
        _, seq_dict, present_taxa, og_len = CreateConcatenationMatrix._process_alignment_file(
            str(fasta),
            ["A", "B", "C"],
        )
        assert present_taxa == {"A", "B"}
        assert og_len == 4
        assert seq_dict["A"] == "ACGT"
        assert seq_dict["B"] == "CCCC"
        assert seq_dict["C"] == "????"

    def test_process_alignment_file_all_present_returns_parsed_sequences(self, tmp_path):
        fasta = tmp_path / "gene.fa"
        _write_fasta(fasta, [("A", "ACGT"), ("B", "CCCC")])

        _, seq_dict, present_taxa, og_len = CreateConcatenationMatrix._process_alignment_file(
            str(fasta),
            ["A", "B"],
        )

        assert present_taxa == {"A", "B"}
        assert og_len == 4
        assert seq_dict == {"A": "ACGT", "B": "CCCC"}

    def test_process_alignment_file_uses_first_header_token(self, tmp_path):
        fasta = tmp_path / "gene.fa"
        fasta.write_text(
            ">A description here\nAC GT\nAA\n"
            ">B another description\nCC\nCC\n"
        )
        _, seq_dict, present_taxa, og_len = CreateConcatenationMatrix._process_alignment_file(
            str(fasta),
            ["A", "B", "C"],
        )
        assert present_taxa == {"A", "B"}
        assert og_len == 6
        assert seq_dict["A"] == "ACGTAA"
        assert seq_dict["B"] == "CCCC"
        assert seq_dict["C"] == "??????"

    def test_create_concatenation_matrix_sequential(self, tmp_path):
        gene1 = tmp_path / "g1.fa"
        gene2 = tmp_path / "g2.fa"
        _write_fasta(gene1, [("A", "ACGT"), ("B", "A-GT")])
        _write_fasta(gene2, [("A", "TT"), ("C", "TA")])

        alignment_list = tmp_path / "alignments.txt"
        alignment_list.write_text(f"{gene1}\n{gene2}\n")
        prefix = tmp_path / "nested" / "concat"

        ccm = CreateConcatenationMatrix(
            Namespace(
                alignment_list=str(alignment_list),
                prefix=str(prefix),
                json=False,
                plot_occupancy=False,
            )
        )
        ccm.create_concatenation_matrix(str(alignment_list), str(prefix))

        fasta_text = Path(f"{prefix}.fa").read_text()
        occupancy_text = Path(f"{prefix}.occupancy").read_text()
        partition_text = Path(f"{prefix}.partition").read_text()

        assert ">A\nACGTTT\n" in fasta_text
        assert ">B\nA-GT??\n" in fasta_text
        assert ">C\n????TA\n" in fasta_text
        assert "AUTO, " in partition_text
        assert str(gene1) in occupancy_text
        assert str(gene2) in occupancy_text

    def test_create_concatenation_matrix_parallel_fallback(self, tmp_path, monkeypatch):
        class FailingExecutor:
            def __init__(self, *_, **__):
                pass

            def __enter__(self):
                raise RuntimeError("disabled")

            def __exit__(self, exc_type, exc, tb):
                return False

        gene_files = []
        for idx, entries in enumerate(
            [
                [("A", "AA"), ("B", "CC")],
                [("A", "GG"), ("C", "TT")],
                [("B", "TA"), ("C", "AT")],
            ]
        ):
            gene = tmp_path / f"g{idx}.fa"
            _write_fasta(gene, entries)
            gene_files.append(gene)

        alignment_list = tmp_path / "alignments_parallel.txt"
        alignment_list.write_text("\n".join(str(p) for p in gene_files) + "\n")
        prefix = tmp_path / "concat_parallel"

        monkeypatch.setattr(ccm_module, "ProcessPoolExecutor", FailingExecutor)
        ccm = CreateConcatenationMatrix(
            Namespace(alignment_list=str(alignment_list), prefix=str(prefix), json=False, plot_occupancy=False)
        )
        ccm.create_concatenation_matrix(str(alignment_list), str(prefix))

        assert Path(f"{prefix}.fa").exists()
        assert Path(f"{prefix}.occupancy").exists()
        assert Path(f"{prefix}.partition").exists()

    def test_sequence_occupancy_states_ascii_and_unicode(self):
        ascii_states = CreateConcatenationMatrix._sequence_occupancy_states(
            "ACGT-?*XxNn"
        )
        unicode_states = CreateConcatenationMatrix._sequence_occupancy_states(
            "A\u03a9-N"
        )

        np.testing.assert_array_equal(
            ascii_states,
            np.array([2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1], dtype=np.uint8),
        )
        np.testing.assert_array_equal(
            unicode_states,
            np.array([2, 2, 1, 1], dtype=np.uint8),
        )

    def test_build_occupancy_state_matrix(self):
        state_matrix, gene_boundaries = (
            CreateConcatenationMatrix._build_occupancy_state_matrix(
                taxa=["A", "B"],
                concatenated_seqs={"A": ["AC", "N?"], "B": ["--", "TT"]},
                present_taxa_by_gene=[{"A", "B"}, {"A"}],
                gene_lengths=[2, 2],
            )
        )

        np.testing.assert_array_equal(
            state_matrix,
            np.array(
                [
                    [2, 2, 1, 1],
                    [1, 1, 0, 0],
                ],
                dtype=np.uint8,
            ),
        )
        assert gene_boundaries == [2, 4]

    def test_build_occupancy_state_matrix_uses_batched_ascii_path(self, monkeypatch):
        def fail_sequence_path(_sequence):
            raise AssertionError("per-sequence path should not be used")

        monkeypatch.setattr(
            CreateConcatenationMatrix,
            "_sequence_occupancy_states",
            staticmethod(fail_sequence_path),
        )

        state_matrix, gene_boundaries = (
            CreateConcatenationMatrix._build_occupancy_state_matrix(
                taxa=["A", "B"],
                concatenated_seqs={"A": ["AC", "N?"], "B": ["--", "TT"]},
                present_taxa_by_gene=[{"A", "B"}, {"A", "B"}],
                gene_lengths=[2, 2],
            )
        )

        np.testing.assert_array_equal(
            state_matrix,
            np.array(
                [
                    [2, 2, 1, 1],
                    [1, 1, 2, 2],
                ],
                dtype=np.uint8,
            ),
        )
        assert gene_boundaries == [2, 4]

    def test_build_occupancy_state_matrix_ignores_unrequested_present_taxa(self):
        state_matrix, gene_boundaries = (
            CreateConcatenationMatrix._build_occupancy_state_matrix(
                taxa=["A"],
                concatenated_seqs={"A": ["AΩ"]},
                present_taxa_by_gene=[{"A", "Z"}],
                gene_lengths=[2],
            )
        )

        np.testing.assert_array_equal(
            state_matrix,
            np.array([[2, 2]], dtype=np.uint8),
        )
        assert gene_boundaries == [2]

    def test_plot_concatenation_occupancy(self, tmp_path, args):
        pytest.importorskip("matplotlib")
        ccm = CreateConcatenationMatrix(args)
        output_file = tmp_path / "occ.png"
        ccm._plot_concatenation_occupancy(
            taxa=["A", "B"],
            alignment_paths=["g1.fa", "g2.fa"],
            concatenated_seqs={"A": ["AC", "GT"], "B": ["A-", "G?"]},
            present_taxa_by_gene=[{"A", "B"}, {"A", "B"}],
            gene_lengths=[2, 2],
            output_file=str(output_file),
        )
        assert output_file.exists()

    def test_plot_concatenation_occupancy_batches_gene_boundaries(
        self, tmp_path, args, monkeypatch
    ):
        pytest.importorskip("matplotlib")
        import matplotlib.axes
        from matplotlib.collections import LineCollection

        ccm = CreateConcatenationMatrix(args)
        output_file = tmp_path / "occ_boundaries.png"
        line_collections = []
        original_add_collection = matplotlib.axes.Axes.add_collection

        def tracking_add_collection(self, collection, *call_args, **kwargs):
            if isinstance(collection, LineCollection):
                line_collections.append(collection)
            return original_add_collection(self, collection, *call_args, **kwargs)

        def fail_axvline(*call_args, **kwargs):
            raise AssertionError("gene boundaries should be batched")

        monkeypatch.setattr(
            matplotlib.axes.Axes,
            "add_collection",
            tracking_add_collection,
        )
        monkeypatch.setattr(matplotlib.axes.Axes, "axvline", fail_axvline)

        ccm._plot_concatenation_occupancy(
            taxa=["A", "B"],
            alignment_paths=["g1.fa", "g2.fa", "g3.fa", "g4.fa"],
            concatenated_seqs={
                "A": ["AC", "GT", "AA", "CC"],
                "B": ["A-", "G?", "TT", "NN"],
            },
            present_taxa_by_gene=[
                {"A", "B"},
                {"A", "B"},
                {"A", "B"},
                {"A", "B"},
            ],
            gene_lengths=[2, 2, 2, 2],
            output_file=str(output_file),
        )

        assert output_file.exists()
        assert len(line_collections) == 1
        assert len(line_collections[0].get_segments()) == 3

    def test_plot_concatenation_occupancy_counts_rows_with_count_nonzero(
        self, tmp_path, args, monkeypatch
    ):
        pytest.importorskip("matplotlib")
        import numpy as real_np

        calls = []
        ccm = CreateConcatenationMatrix(args)

        def tracking_count_nonzero(values, *call_args, **kwargs):
            calls.append(kwargs.get("axis"))
            return real_np.count_nonzero(values, *call_args, **kwargs)

        def fail_sum(*_args, **_kwargs):
            raise AssertionError("represented occupancy rows should use count_nonzero")

        monkeypatch.setattr(ccm_module.np, "count_nonzero", tracking_count_nonzero)
        monkeypatch.setattr(ccm_module.np, "sum", fail_sum, raising=False)

        ccm._plot_concatenation_occupancy(
            taxa=["A", "B"],
            alignment_paths=["g1.fa", "g2.fa"],
            concatenated_seqs={"A": ["AC", "GT"], "B": ["A-", "G?"]},
            present_taxa_by_gene=[{"A", "B"}, {"A", "B"}],
            gene_lengths=[2, 2],
            output_file=str(tmp_path / "occ_count_nonzero.png"),
        )

        assert calls == [1]

    def test_plot_concatenation_occupancy_pdf_output(self, tmp_path, args):
        pytest.importorskip("matplotlib")
        ccm = CreateConcatenationMatrix(args)
        output_file = tmp_path / "occ.pdf"
        ccm._plot_concatenation_occupancy(
            taxa=["A", "B"],
            alignment_paths=["g1.fa", "g2.fa"],
            concatenated_seqs={"A": ["AC", "GT"], "B": ["A-", "G?"]},
            present_taxa_by_gene=[{"A", "B"}, {"A", "B"}],
            gene_lengths=[2, 2],
            output_file=str(output_file),
        )
        assert output_file.exists()

    def test_plot_concatenation_occupancy_custom_config(self, tmp_path):
        pytest.importorskip("matplotlib")
        custom_args = Namespace(
            alignment_list="/x", prefix="x",
            fig_width=10.0, fig_height=6.0, dpi=72, no_title=True,
            title=None, legend_position="lower left",
            ylabel_fontsize=5.0, xlabel_fontsize=4.0,
            title_fontsize=10.0, axis_fontsize=8.0,
            colors="#000000,#ffffff,#ff0000",
        )
        ccm = CreateConcatenationMatrix(custom_args)
        output_file = tmp_path / "custom.png"
        ccm._plot_concatenation_occupancy(
            taxa=["A", "B", "C"],
            alignment_paths=["g1.fa", "g2.fa"],
            concatenated_seqs={"A": ["AC", "GT"], "B": ["A-", "G?"], "C": ["TT", "AA"]},
            present_taxa_by_gene=[{"A", "B"}, {"A", "C"}],
            gene_lengths=[2, 2],
            output_file=str(output_file),
        )
        assert output_file.exists()

    def test_plot_concatenation_occupancy_svg_output(self, tmp_path, args):
        pytest.importorskip("matplotlib")
        ccm = CreateConcatenationMatrix(args)
        output_file = tmp_path / "occ.svg"
        ccm._plot_concatenation_occupancy(
            taxa=["A", "B"],
            alignment_paths=["g1.fa", "g2.fa"],
            concatenated_seqs={"A": ["AC", "GT"], "B": ["A-", "G?"]},
            present_taxa_by_gene=[{"A", "B"}, {"A", "B"}],
            gene_lengths=[2, 2],
            output_file=str(output_file),
        )
        assert output_file.exists()

    def test_create_concatenation_matrix_json_and_plot(self, tmp_path, monkeypatch):
        gene1 = tmp_path / "g1.fa"
        gene2 = tmp_path / "g2.fa"
        _write_fasta(gene1, [("A", "AC"), ("B", "TT")])
        _write_fasta(gene2, [("A", "GA"), ("C", "TA")])
        alignment_list = tmp_path / "alignments.txt"
        alignment_list.write_text(f"{gene1}\n{gene2}\n")
        prefix = tmp_path / "concat_json"
        captured = {}

        ccm = CreateConcatenationMatrix(
            Namespace(
                alignment_list=str(alignment_list),
                prefix=str(prefix),
                json=True,
                plot_occupancy=True,
                plot_output=str(tmp_path / "custom_plot.png"),
            )
        )

        monkeypatch.setattr(ccm, "_plot_concatenation_occupancy", lambda **kwargs: captured.setdefault("plot", kwargs))
        monkeypatch.setattr(ccm_module, "print_json", lambda payload: captured.setdefault("payload", payload))
        ccm.create_concatenation_matrix(str(alignment_list), str(prefix))

        assert captured["plot"]["output_file"] == str(tmp_path / "custom_plot.png")
        assert captured["payload"]["total_alignments"] == 2
        assert captured["payload"]["total_taxa"] == 3
        assert captured["payload"]["concatenated_length"] == 4
        assert "occupancy_plot" in captured["payload"]["output_files"]

    def test_compute_effective_occupancy_all_valid(self, args):
        ccm = CreateConcatenationMatrix(args)
        seqs = {"A": ["ACGT", "TTGG"]}
        occupancy, excluded = ccm._compute_effective_occupancy(seqs, 0.5)
        assert occupancy["A"] == 1.0
        assert excluded == set()

    def test_compute_effective_occupancy_all_gaps(self, args):
        ccm = CreateConcatenationMatrix(args)
        seqs = {"A": ["----", "NNNN"]}
        occupancy, excluded = ccm._compute_effective_occupancy(seqs, 0.5)
        assert occupancy["A"] == 0.0
        assert "A" in excluded

    def test_compute_effective_occupancy_mixed(self, args):
        ccm = CreateConcatenationMatrix(args)
        # 4 informative out of 8 total = 0.5
        seqs = {"A": ["AC--", "??GT"]}
        occupancy, excluded = ccm._compute_effective_occupancy(seqs, 0.5)
        assert occupancy["A"] == 0.5
        # 0.5 is not < 0.5, so not excluded
        assert excluded == set()

    def test_compute_effective_occupancy_counts_all_invalid_symbols(self, args):
        ccm = CreateConcatenationMatrix(args)
        seqs = {"A": ["ACGT-?*XxNn", "AA"]}
        occupancy, excluded = ccm._compute_effective_occupancy(seqs, 0.5)
        assert occupancy["A"] == pytest.approx(6 / 13)
        assert excluded == {"A"}

    def test_compute_effective_occupancy_unicode_fallback(self, args):
        ccm = CreateConcatenationMatrix(args)
        seqs = {"A": ["A\u03a9-N", "??"]}

        occupancy, excluded = ccm._compute_effective_occupancy(seqs, 0.5)

        assert occupancy["A"] == pytest.approx(2 / 6)
        assert excluded == {"A"}

    def test_compute_effective_occupancy_counts_across_many_parts(self, args):
        ccm = CreateConcatenationMatrix(args)
        seqs = {
            "A": ["AC", "", "GT", "NN"],
            "B": ["--", "??", "*X", "xn"],
        }

        occupancy, excluded = ccm._compute_effective_occupancy(seqs, 0.5)

        assert occupancy["A"] == pytest.approx(4 / 6)
        assert occupancy["B"] == 0.0
        assert excluded == {"B"}

    def test_threshold_excludes_low_occupancy_taxon(self, tmp_path):
        gene1 = tmp_path / "g1.fa"
        gene2 = tmp_path / "g2.fa"
        # A is fully present in both genes
        # B has only gaps/ambiguous chars
        _write_fasta(gene1, [("A", "ACGT"), ("B", "----")])
        _write_fasta(gene2, [("A", "TTGG"), ("B", "NNNN")])

        alignment_list = tmp_path / "alignments.txt"
        alignment_list.write_text(f"{gene1}\n{gene2}\n")
        prefix = tmp_path / "concat_thresh"

        ccm = CreateConcatenationMatrix(
            Namespace(
                alignment_list=str(alignment_list),
                prefix=str(prefix),
                json=False,
                plot_occupancy=False,
                threshold=0.5,
            )
        )
        ccm.create_concatenation_matrix(str(alignment_list), str(prefix))

        fasta_text = Path(f"{prefix}.fa").read_text()
        assert ">A\n" in fasta_text
        assert ">B\n" not in fasta_text

    def test_threshold_report_batches_excluded_taxa(self, tmp_path, mocker):
        ccm = CreateConcatenationMatrix(
            Namespace(
                alignment_list="alignments.txt",
                prefix=str(tmp_path / "concat"),
                json=False,
                plot_occupancy=False,
                threshold=0.5,
            )
        )
        records = [
            _ParsedFastaRecord("A", "ACGT"),
            _ParsedFastaRecord("B", "----"),
            _ParsedFastaRecord("C", "NNNN"),
        ]

        mocker.patch.object(ccm, "read_alignment_paths", return_value=["g1.fa"])
        mocker.patch.object(ccm, "get_taxa_names", return_value=["A", "B", "C"])
        mocker.patch.object(
            ccm,
            "get_list_of_taxa_and_records",
            return_value=({"A", "B", "C"}, records),
        )
        mocker.patch.object(ccm, "print_start_message")
        mocker.patch.object(ccm, "fasta_file_write")
        mocker.patch.object(ccm, "write_occupancy_or_partition_file")
        printed = mocker.patch("builtins.print")

        ccm.create_concatenation_matrix("alignments.txt", str(tmp_path / "concat"))

        assert printed.call_args_list[0].args == (
            "Threshold (0.5): excluded 2 taxa\n"
            "  B: effective occupancy = 0.0\n"
            "  C: effective occupancy = 0.0",
        )
        assert printed.call_args_list[1].args == ("Complete!\n",)
        assert printed.call_count == 2

    def test_threshold_zero_disables_filtering(self, tmp_path):
        gene1 = tmp_path / "g1.fa"
        gene2 = tmp_path / "g2.fa"
        _write_fasta(gene1, [("A", "ACGT"), ("B", "----")])
        _write_fasta(gene2, [("A", "TTGG"), ("B", "NNNN")])

        alignment_list = tmp_path / "alignments.txt"
        alignment_list.write_text(f"{gene1}\n{gene2}\n")
        prefix = tmp_path / "concat_nofilt"

        ccm = CreateConcatenationMatrix(
            Namespace(
                alignment_list=str(alignment_list),
                prefix=str(prefix),
                json=False,
                plot_occupancy=False,
                threshold=0,
            )
        )
        ccm.create_concatenation_matrix(str(alignment_list), str(prefix))

        fasta_text = Path(f"{prefix}.fa").read_text()
        assert ">A\n" in fasta_text
        assert ">B\n" in fasta_text

    def test_threshold_one_excludes_all_but_complete(self, tmp_path):
        gene1 = tmp_path / "g1.fa"
        gene2 = tmp_path / "g2.fa"
        # A: fully present in both → occupancy 1.0
        # B: present in gene1, missing from gene2 → has ???? → occupancy 0.5
        _write_fasta(gene1, [("A", "ACGT"), ("B", "TTGG")])
        _write_fasta(gene2, [("A", "AACC")])

        alignment_list = tmp_path / "alignments.txt"
        alignment_list.write_text(f"{gene1}\n{gene2}\n")
        prefix = tmp_path / "concat_strict"

        ccm = CreateConcatenationMatrix(
            Namespace(
                alignment_list=str(alignment_list),
                prefix=str(prefix),
                json=False,
                plot_occupancy=False,
                threshold=1.0,
            )
        )
        ccm.create_concatenation_matrix(str(alignment_list), str(prefix))

        fasta_text = Path(f"{prefix}.fa").read_text()
        assert ">A\n" in fasta_text
        assert ">B\n" not in fasta_text
