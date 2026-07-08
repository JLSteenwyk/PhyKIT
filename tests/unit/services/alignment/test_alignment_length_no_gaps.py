import pytest
import subprocess
import sys
from argparse import Namespace
from math import isclose
from types import SimpleNamespace

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.alignment_length_no_gaps import AlignmentLengthNoGaps
import phykit.services.alignment.alignment_length_no_gaps as alg_module


@pytest.fixture
def args():
    kwargs = dict(alignment="/some/path/to/file.fa")
    return Namespace(**kwargs)


class TestAlignmentLengthNoGaps(object):
    def test_module_import_does_not_import_numpy_or_biopython_align(self):
        code = """
import sys
import phykit.services.alignment.alignment_length_no_gaps as module
assert hasattr(module.np, "__getattr__")
assert "argparse" not in sys.modules
assert "typing" not in sys.modules
assert "numpy" not in sys.modules
assert "Bio.Align" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
"""
        subprocess.run([sys.executable, "-c", code], check=True)

    def test_lazy_numpy_caches_resolved_attributes(self):
        lazy_np = alg_module._LazyNumpy()

        any_attr = lazy_np.any

        assert lazy_np.__dict__["any"] is any_attr
        assert lazy_np.any is any_attr
        assert lazy_np._module is not None

    def test_all_sequences_identical_does_not_slice(self):
        class NoSliceList(list):
            def __getitem__(self, key):
                if isinstance(key, slice):
                    raise AssertionError("identical sequence scan should not slice")
                return super().__getitem__(key)

        assert alg_module._all_sequences_identical(
            NoSliceList(["ACGT", "ACGT", "ACGT"])
        )
        assert not alg_module._all_sequences_identical(
            NoSliceList(["ACGT", "TGCA", "ACGT"])
        )
        assert not alg_module._all_sequences_identical(
            NoSliceList(["ACGT", "ACGT", "TGCA"])
        )

    def test_init_sets_alignment_file_path(self, args):
        aln = AlignmentLengthNoGaps(args)
        assert aln.alignment_file_path == args.alignment
        assert aln.output_file_path is None

    def test_alignment_length_no_gaps(self, alignment_simple, args):
        aln = AlignmentLengthNoGaps(args)
        is_protein = True
        (
            aln_len_no_gaps,
            aln_len,
            aln_len_no_gaps_per,
        ) = aln.calculate_alignment_length_no_gaps(alignment_simple, is_protein)
        assert isinstance(aln_len_no_gaps, int)
        assert isinstance(aln_len, int)
        assert isinstance(aln_len_no_gaps_per, float)
        assert aln_len_no_gaps == 3
        assert aln_len == 6
        assert isclose(aln_len_no_gaps_per, 50.0, rel_tol=0.001)


    def test_calculate_alignment_length_no_gaps(self, mocker, args):
        expected_length = 6
        expected_aln_len_no_gaps = 3
        expected_aln_len_no_gaps_per = 50.0
        mocker.patch(
            "phykit.services.alignment.alignment_length_no_gaps.AlignmentLengthNoGaps.get_sites_no_gaps_count",
            return_value=expected_aln_len_no_gaps
        )
        mock_aln = mocker.MagicMock(
            get_alignment_length=mocker.MagicMock(return_value=expected_length)
        )

        aln = AlignmentLengthNoGaps(args) 
        res = aln.calculate_alignment_length_no_gaps(mock_aln, True)
        assert res[0] == expected_aln_len_no_gaps
        assert res[1] == expected_length
        assert res[2] == expected_aln_len_no_gaps_per

    def test_get_sites_no_gaps_count(self, mocker, args, alignment_complex):
        aln = AlignmentLengthNoGaps(args)
        is_protein = True
        expected_length = 955
        expected_aln_len_no_gaps = 901
        res = aln.get_sites_no_gaps_count(alignment_complex, expected_length, is_protein)
        assert res == expected_aln_len_no_gaps

    def test_get_sites_no_gaps_count_handles_lowercase_ambiguity(self, args):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AcGtA"), id="a"),
                SeqRecord(Seq("AnG-A"), id="b"),
                SeqRecord(Seq("aCGtx"), id="c"),
            ]
        )
        aln = AlignmentLengthNoGaps(args)

        dna_count = aln.get_sites_no_gaps_count(alignment, 5, is_protein=False)
        protein_count = aln.get_sites_no_gaps_count(alignment, 5, is_protein=True)

        assert dna_count == 2
        assert protein_count == 3

    def test_get_sites_no_gaps_count_ascii_path_avoids_isin(self, mocker, args):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AcGtA"), id="a"),
                SeqRecord(Seq("AnG-A"), id="b"),
                SeqRecord(Seq("aCGtx"), id="c"),
            ]
        )
        aln = AlignmentLengthNoGaps(args)
        mocker.patch(
            "phykit.services.alignment.alignment_length_no_gaps.np.isin",
            side_effect=AssertionError("ASCII path should use the lookup table"),
        )

        assert aln.get_sites_no_gaps_count(alignment, 5, is_protein=False) == 2

    def test_get_sites_no_gaps_count_ascii_path_uses_gap_byte_positions(
        self, mocker, args
    ):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGTA"), id="a"),
                SeqRecord(Seq("ANG-A"), id="b"),
                SeqRecord(Seq("ACGTx"), id="c"),
            ]
        )
        aln = AlignmentLengthNoGaps(args)
        gap_column_spy = mocker.spy(alg_module, "_count_columns_without_gap_bytes")
        mocker.patch(
            "phykit.services.alignment.alignment_length_no_gaps.np.frombuffer",
            side_effect=AssertionError("ASCII gap path should scan gap byte positions"),
        )

        assert aln.get_sites_no_gaps_count(alignment, 5, is_protein=False) == 2
        gap_column_spy.assert_called_once()

    def test_count_columns_without_gap_bytes_short_circuits_full_gap_coverage(self):
        alignment_bytes = b"-N-N-"

        assert (
            alg_module._count_columns_without_gap_bytes(
                alignment_bytes,
                aln_len=5,
                gap_bytes=b"-N",
            )
            == 0
        )

    def test_get_sites_no_gaps_count_ascii_path_skips_gap_char_set(
        self, mocker, args
    ):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGTA"), id="a"),
                SeqRecord(Seq("ANG-A"), id="b"),
                SeqRecord(Seq("ACGTx"), id="c"),
            ]
        )
        aln = AlignmentLengthNoGaps(args)
        mocker.patch.object(
            aln,
            "get_gap_chars",
            side_effect=AssertionError(
                "ASCII path should not build Unicode fallback gap characters"
            ),
        )

        assert aln.get_sites_no_gaps_count(alignment, 5, is_protein=False) == 2

    def test_get_sites_no_gaps_count_no_gap_ascii_returns_length(self, mocker, args):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGTA"), id="a"),
                SeqRecord(Seq("ACGTT"), id="b"),
                SeqRecord(Seq("ACGTC"), id="c"),
            ]
        )
        aln = AlignmentLengthNoGaps(args)
        mocker.patch(
            "phykit.services.alignment.alignment_length_no_gaps.np.any",
            side_effect=AssertionError("no-gap ASCII path should not reduce columns"),
        )
        mocker.patch(
            "phykit.services.alignment.alignment_length_no_gaps."
            "_count_columns_without_gap_bytes",
            side_effect=AssertionError("no-gap ASCII path should not scan gap columns"),
        )

        assert aln.get_sites_no_gaps_count(alignment, 5, is_protein=False) == 5

    def test_get_sites_no_gaps_count_identical_sequences_skip_matrix(
        self, mocker, args
    ):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGTNn-?*Xx"), id="a"),
                SeqRecord(Seq("ACGTNn-?*Xx"), id="b"),
                SeqRecord(Seq("ACGTNn-?*Xx"), id="c"),
            ]
        )
        aln = AlignmentLengthNoGaps(args)
        mocker.patch(
            "phykit.services.alignment.alignment_length_no_gaps.np.frombuffer",
            side_effect=AssertionError(
                "Identical sequences should avoid alignment matrix construction"
            ),
        )

        assert aln.get_sites_no_gaps_count(alignment, 11, is_protein=False) == 4
        assert aln.get_sites_no_gaps_count(alignment, 11, is_protein=True) == 6

    def test_get_sites_no_gaps_count_identical_unicode_sequence(self, args):
        alignment = [
            SimpleNamespace(seq="A\u03a9n-?", id="a"),
            SimpleNamespace(seq="A\u03a9n-?", id="b"),
        ]
        aln = AlignmentLengthNoGaps(args)

        assert aln.get_sites_no_gaps_count(alignment, 5, is_protein=False) == 2
        assert aln.get_sites_no_gaps_count(alignment, 5, is_protein=True) == 3

    def test_get_sites_no_gaps_count_identical_unicode_lowercase_gap_codes(
        self, args
    ):
        alignment = [
            SimpleNamespace(seq="A\u03a9xN-?", id="a"),
            SimpleNamespace(seq="A\u03a9xN-?", id="b"),
        ]
        aln = AlignmentLengthNoGaps(args)

        assert aln.get_sites_no_gaps_count(alignment, 6, is_protein=False) == 2
        assert aln.get_sites_no_gaps_count(alignment, 6, is_protein=True) == 3

    def test_get_sites_no_gaps_count_unicode_fallback(self, args):
        alignment = [
            SimpleNamespace(seq="A\u03a9nT", id="a"),
            SimpleNamespace(seq="ACGT", id="b"),
            SimpleNamespace(seq="ATGT", id="c"),
        ]
        aln = AlignmentLengthNoGaps(args)

        assert aln.get_sites_no_gaps_count(alignment, 4, is_protein=False) == 3
        assert aln.get_sites_no_gaps_count(alignment, 4, is_protein=True) == 4
