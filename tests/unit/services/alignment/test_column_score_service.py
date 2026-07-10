from argparse import Namespace
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pytest
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.column_score import ColumnScore
import phykit.services.alignment.column_score as module


def test_module_import_does_not_import_numpy_bio_align_or_json_helpers():
    code = """
import sys
import phykit.services.alignment.column_score as module
assert callable(module.get_alignment_and_format)
assert hasattr(module.AlignIO, "read")
assert hasattr(module.np, "__getattr__")
assert "typing" not in sys.modules
assert "numpy" not in sys.modules
assert "Bio.Align" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_lazy_numpy_caches_resolved_attributes():
    lazy_np = module._LazyNumpy()

    frombuffer_attr = lazy_np.frombuffer

    assert lazy_np.__dict__["frombuffer"] is frombuffer_attr
    assert lazy_np.frombuffer is frombuffer_attr
    assert lazy_np._module is not None


def test_lazy_alignio_caches_resolved_read():
    lazy_alignio = module._LazyAlignIO()

    from Bio import AlignIO as real_alignio

    original_read = real_alignio.read
    calls = []

    def fake_read(*args, **kwargs):
        calls.append((args, kwargs))
        return "alignment"

    try:
        real_alignio.read = fake_read

        assert lazy_alignio.read("query.fa", "fasta") == "alignment"
        cached_read = lazy_alignio.__dict__["read"]

        def changed_read(*args, **kwargs):
            return "changed"

        real_alignio.read = changed_read

        assert lazy_alignio.read("reference.fa", "fasta") == "alignment"
        assert lazy_alignio.__dict__["read"] is cached_read
        assert calls == [
            (("query.fa", "fasta"), {}),
            (("reference.fa", "fasta"), {}),
        ]
    finally:
        real_alignio.read = original_read


@pytest.fixture
def args():
    return Namespace(fasta="/some/path/to/query.fa", reference="/some/path/to/ref.fa")


class TestColumnScore:
    def test_init_sets_expected_attrs(self, args):
        service = ColumnScore(args)
        assert service.fasta == args.fasta
        assert service.reference == args.reference
        assert service.json_output is False

    def test_get_columns_from_alignments(self, args):
        service = ColumnScore(args)
        reference = MultipleSeqAlignment(
            [SeqRecord(Seq("aC"), id="r1"), SeqRecord(Seq("Gt"), id="r2")]
        )
        query = MultipleSeqAlignment(
            [SeqRecord(Seq("At"), id="q1"), SeqRecord(Seq("gT"), id="q2")]
        )
        ref_cols, query_cols = service.get_columns_from_alignments(reference, query)
        assert ref_cols == ["AG", "CT"]
        assert query_cols == ["AG", "TT"]

    def test_calculate_matches_between_ref_and_query_columns(self, args):
        service = ColumnScore(args)
        matches, total = service.calculate_matches_between_ref_and_query_columns(
            ["AA", "CC", "GG"],
            ["AA", "TT", "GG"],
        )
        assert matches == 2
        assert total == 3

    def test_direct_column_matching_matches_string_helper_for_ascii(self, args):
        service = ColumnScore(args)
        reference = MultipleSeqAlignment(
            [SeqRecord(Seq("aCT"), id="r1"), SeqRecord(Seq("GtT"), id="r2")]
        )
        query = MultipleSeqAlignment(
            [SeqRecord(Seq("AtA"), id="q1"), SeqRecord(Seq("gTT"), id="q2")]
        )

        ref_cols, query_cols = service.get_columns_from_alignments(reference, query)
        expected = service.calculate_matches_between_ref_and_query_columns(
            ref_cols, query_cols
        )

        assert service._calculate_matches_between_alignments_direct(
            reference, query
        ) == expected

    def test_direct_column_matching_same_object_uniques_once_preserving_duplicates(
        self, mocker, args
    ):
        service = ColumnScore(args)
        alignment = MultipleSeqAlignment(
            [SeqRecord(Seq("AAAA"), id="r1"), SeqRecord(Seq("CCCC"), id="r2")]
        )
        intersect_spy = mocker.spy(np, "intersect1d")

        assert service._calculate_matches_between_alignments_direct(
            alignment, alignment
        ) == (1, 4)
        intersect_spy.assert_not_called()

    def test_direct_column_matching_same_object_repeated_rows_skips_matrix(
        self, mocker, args
    ):
        service = ColumnScore(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGTNA"), id="r1"),
                SeqRecord(Seq("acgtna"), id="r2"),
                SeqRecord(Seq("ACGTNA"), id="r3"),
            ]
        )
        mocker.patch(
            "phykit.services.alignment.column_score.np.frombuffer",
            side_effect=AssertionError("repeated rows should not build a matrix"),
        )

        assert service._calculate_matches_between_alignments_direct(
            alignment, alignment
        ) == (5, 6)

    def test_direct_column_matching_identical_sequences_uniques_once(
        self, mocker, args
    ):
        service = ColumnScore(args)
        reference = MultipleSeqAlignment(
            [SeqRecord(Seq("AAAA"), id="r1"), SeqRecord(Seq("CCCC"), id="r2")]
        )
        query = MultipleSeqAlignment(
            [SeqRecord(Seq("aaaa"), id="q1"), SeqRecord(Seq("cccc"), id="q2")]
        )
        intersect_spy = mocker.spy(np, "intersect1d")

        assert service._calculate_matches_between_alignments_direct(
            reference, query
        ) == (1, 4)
        intersect_spy.assert_not_called()

    def test_direct_column_matching_small_taxon_sets_skip_numpy_matrix(
        self, mocker, args
    ):
        service = ColumnScore(args)
        reference = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGTAC"), id="r1"),
                SeqRecord(Seq("TGCATG"), id="r2"),
                SeqRecord(Seq("AAAAAA"), id="r3"),
                SeqRecord(Seq("CCCCCC"), id="r4"),
            ]
        )
        query = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGTAA"), id="q1"),
                SeqRecord(Seq("TGCATT"), id="q2"),
                SeqRecord(Seq("AAAACC"), id="q3"),
                SeqRecord(Seq("CCCCGG"), id="q4"),
            ]
        )
        mocker.patch(
            "phykit.services.alignment.column_score.np.frombuffer",
            side_effect=AssertionError("small taxon sets should use column sets"),
        )

        assert service._calculate_matches_between_alignments_direct(
            reference, query
        ) == (4, 6)

    def test_direct_column_matching_five_taxa_skip_numpy_matrix(
        self, mocker, args
    ):
        service = ColumnScore(args)
        reference = MultipleSeqAlignment(
            [
                SeqRecord(Seq("A-G-TAT"), id="r1"),
                SeqRecord(Seq("A-G--AT"), id="r2"),
                SeqRecord(Seq("A-G--TA"), id="r3"),
                SeqRecord(Seq("AG-A-TA"), id="r4"),
                SeqRecord(Seq("AC-a-T-"), id="r5"),
            ]
        )
        query = MultipleSeqAlignment(
            [
                SeqRecord(Seq("A-GTAT"), id="q1"),
                SeqRecord(Seq("A-G-AT"), id="q2"),
                SeqRecord(Seq("A-G-TA"), id="q3"),
                SeqRecord(Seq("AGA-TA"), id="q4"),
                SeqRecord(Seq("ACa-T-"), id="q5"),
            ]
        )
        mocker.patch(
            "phykit.services.alignment.column_score.np.frombuffer",
            side_effect=AssertionError("five-taxon column score should use sets"),
        )

        assert service._calculate_matches_between_alignments_direct(
            reference, query
        ) == (5, 6)

    def test_small_ascii_column_set_rejects_unicode(self):
        assert ColumnScore._small_ascii_column_set(["AΩ", "AT"]) is None

    def test_small_ascii_column_set_single_taxon_uses_symbol_set(self):
        assert ColumnScore._small_ascii_column_set(["ACGTA"]) == {"A", "C", "G", "T"}

    def test_direct_column_matching_single_taxon_intersects_symbols(
        self, mocker, args
    ):
        service = ColumnScore(args)
        reference = MultipleSeqAlignment([SeqRecord(Seq("ACGTAC"), id="r1")])
        query = MultipleSeqAlignment([SeqRecord(Seq("CCTTAA"), id="q1")])
        mocker.patch(
            "phykit.services.alignment.column_score.np.frombuffer",
            side_effect=AssertionError("single-taxon columns should use symbols"),
        )

        assert service._calculate_matches_between_alignments_direct(
            reference, query
        ) == (3, 6)

    def test_direct_column_matching_repeated_rows_intersects_symbols(
        self, mocker, args
    ):
        service = ColumnScore(args)
        reference = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AACCGG"), id="r1"),
                SeqRecord(Seq("AACCGG"), id="r2"),
            ]
        )
        query = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ccttaa"), id="q1"),
                SeqRecord(Seq("CCTTAA"), id="q2"),
            ]
        )
        mocker.patch(
            "phykit.services.alignment.column_score.np.frombuffer",
            side_effect=AssertionError("repeated rows should not build a matrix"),
        )

        assert service._calculate_matches_between_alignments_direct(
            reference, query
        ) == (2, 6)

    def test_repeated_sequence_symbols_ascii_returns_symbols_for_matching_rows(self):
        assert ColumnScore._repeated_sequence_symbols_ascii(
            ["AAGT", "AAGT", "AAGT"]
        ) == frozenset({"A", "G", "T"})

    def test_direct_column_matching_handles_taxon_count_mismatch(self, args):
        service = ColumnScore(args)
        reference = MultipleSeqAlignment(
            [SeqRecord(Seq("AC"), id="r1"), SeqRecord(Seq("GT"), id="r2")]
        )
        query = MultipleSeqAlignment([SeqRecord(Seq("AC"), id="q1")])

        assert service._calculate_matches_between_alignments_direct(
            reference, query
        ) == (0, 2)

    def test_direct_column_matching_taxon_count_mismatch_skips_sequence_materialization(
        self, args
    ):
        class UnstringableSequence:
            def __str__(self):
                raise AssertionError(
                    "taxon-count mismatch should not inspect sequences"
                )

        class DummyAlignment(list):
            def get_alignment_length(self):
                return 5

        service = ColumnScore(args)
        reference = DummyAlignment(
            [
                SimpleNamespace(seq=UnstringableSequence()),
                SimpleNamespace(seq=UnstringableSequence()),
            ]
        )
        query = DummyAlignment([SimpleNamespace(seq=UnstringableSequence())])

        assert service._calculate_matches_between_alignments_direct(
            reference,
            query,
        ) == (0, 5)

    def test_calculate_fasta_paths_direct_uses_sequence_lists(self, tmp_path):
        query = tmp_path / "query.fa"
        reference = tmp_path / "reference.fa"
        query.write_text(
            ">1\nA-GTAT\n"
            ">2\nA-G-AT\n"
            ">3\nA-G-TA\n"
            ">4\nAGA-TA\n"
            ">5\nACa-T-\n"
        )
        reference.write_text(
            ">1\nA-G-TAT\n"
            ">2\nA-G--AT\n"
            ">3\nA-G--TA\n"
            ">4\nAG-A-TA\n"
            ">5\nAC-a-T-\n"
        )

        assert ColumnScore._calculate_fasta_paths_direct(
            str(reference),
            str(query),
        ) == (5, 6)

    def test_calculate_fasta_paths_direct_falls_back_for_non_fasta(self, tmp_path):
        query = tmp_path / "query.phy"
        reference = tmp_path / "reference.fa"
        query.write_text("2 4\nr1 ACGT\nr2 ACGT\n")
        reference.write_text(">1\nACGT\n")

        assert ColumnScore._calculate_fasta_paths_direct(
            str(reference),
            str(query),
        ) is None

    def test_run_uses_direct_fasta_path(self, mocker, tmp_path):
        query = tmp_path / "query.fa"
        reference = tmp_path / "reference.fa"
        query.write_text(">1\nACGTAC\n>2\nTGCATG\n")
        reference.write_text(">1\nACGTAA\n>2\nTGCATT\n")
        args = Namespace(fasta=str(query), reference=str(reference), json=False)
        service = ColumnScore(args)
        mocker.patch(
            "phykit.services.alignment.column_score.get_alignment_and_format",
            side_effect=AssertionError("simple FASTA should use the direct path"),
        )
        mocked_print = mocker.patch("builtins.print")

        service.run()

        mocked_print.assert_called_once_with(0.6667)

    def test_direct_column_matching_falls_back_for_unicode(self, args):
        service = ColumnScore(args)
        reference = [type("Record", (), {"seq": "AΩ"})()]
        query = [type("Record", (), {"seq": "AΩ"})()]

        assert service._calculate_matches_between_alignments_direct(
            reference, query
        ) is None

    def test_run_prints_score(self, mocker):
        args = Namespace(fasta="/some/path/to/query.fa", reference="/some/path/to/ref.fa", json=False)
        service = ColumnScore(args)
        query = MultipleSeqAlignment(
            [SeqRecord(Seq("AT"), id="q1"), SeqRecord(Seq("GT"), id="q2")]
        )
        reference = MultipleSeqAlignment(
            [SeqRecord(Seq("AC"), id="r1"), SeqRecord(Seq("GT"), id="r2")]
        )
        mocker.patch(
            "phykit.services.alignment.column_score.get_alignment_and_format",
            side_effect=[(query, "fasta", False), (reference, "fasta", False)],
        )
        mocked_print = mocker.patch("builtins.print")
        service.run()
        mocked_print.assert_called_once_with(0.5)

    def test_run_uses_direct_column_matching(self, mocker):
        args = Namespace(fasta="/some/path/to/query.fa", reference="/some/path/to/ref.fa", json=False)
        service = ColumnScore(args)
        query = MultipleSeqAlignment(
            [SeqRecord(Seq("AT"), id="q1"), SeqRecord(Seq("GT"), id="q2")]
        )
        reference = MultipleSeqAlignment(
            [SeqRecord(Seq("AC"), id="r1"), SeqRecord(Seq("GT"), id="r2")]
        )
        mocker.patch(
            "phykit.services.alignment.column_score.get_alignment_and_format",
            side_effect=[(query, "fasta", False), (reference, "fasta", False)],
        )
        mocker.patch.object(
            ColumnScore,
            "get_columns_from_alignments",
            side_effect=AssertionError("ASCII alignments should use the direct path"),
        )
        mocked_print = mocker.patch("builtins.print")

        service.run()

        mocked_print.assert_called_once_with(0.5)

    def test_run_reads_once_when_query_and_reference_paths_match(self, mocker):
        args = Namespace(fasta="/same/path.fa", reference="/same/path.fa", json=False)
        service = ColumnScore(args)
        alignment = MultipleSeqAlignment(
            [SeqRecord(Seq("AAAA"), id="r1"), SeqRecord(Seq("CCCC"), id="r2")]
        )
        read_mock = mocker.patch(
            "phykit.services.alignment.column_score.get_alignment_and_format",
            return_value=(alignment, "fasta", False),
        )
        mocked_print = mocker.patch("builtins.print")

        service.run()

        read_mock.assert_called_once_with("/same/path.fa")
        mocked_print.assert_called_once_with(0.25)

    def test_run_falls_back_to_alignio_when_helper_detects_non_fasta(self, mocker):
        args = Namespace(fasta="/some/path/to/query.fa", reference="/some/path/to/ref.fa", json=False)
        service = ColumnScore(args)
        query = MultipleSeqAlignment(
            [SeqRecord(Seq("AT"), id="q1"), SeqRecord(Seq("GT"), id="q2")]
        )
        reference = MultipleSeqAlignment(
            [SeqRecord(Seq("AC"), id="r1"), SeqRecord(Seq("GT"), id="r2")]
        )
        mocker.patch(
            "phykit.services.alignment.column_score.get_alignment_and_format",
            side_effect=[(object(), "phylip", False), (object(), "phylip", False)],
        )
        read_mock = mocker.patch(
            "phykit.services.alignment.column_score.AlignIO.read",
            side_effect=[query, reference],
        )
        mocked_print = mocker.patch("builtins.print")

        service.run()

        assert read_mock.call_args_list == [
            mocker.call("/some/path/to/query.fa", "fasta"),
            mocker.call("/some/path/to/ref.fa", "fasta"),
        ]
        mocked_print.assert_called_once_with(0.5)

    def test_run_json_output(self, mocker):
        args = Namespace(fasta="/some/path/to/query.fa", reference="/some/path/to/ref.fa", json=True)
        service = ColumnScore(args)
        query = MultipleSeqAlignment(
            [SeqRecord(Seq("AT"), id="q1"), SeqRecord(Seq("GT"), id="q2")]
        )
        reference = MultipleSeqAlignment(
            [SeqRecord(Seq("AC"), id="r1"), SeqRecord(Seq("GT"), id="r2")]
        )
        mocker.patch(
            "phykit.services.alignment.column_score.get_alignment_and_format",
            side_effect=[(query, "fasta", False), (reference, "fasta", False)],
        )
        mocked_json = mocker.patch("phykit.services.alignment.column_score.print_json")
        service.run()
        mocked_json.assert_called_once_with({"column_score": 0.5})
