import pytest
import subprocess
import sys
from argparse import Namespace
from math import isclose
from types import SimpleNamespace

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.rcv import RelativeCompositionVariability
import phykit.services.alignment.rcv as rcv_module
import phykit.services.alignment.base as alignment_base_module


@pytest.fixture
def args():
    kwargs = dict(alignment="/some/path/to/file.fa")
    return Namespace(**kwargs)


class TestRelativeCompositionVariability(object):
    def test_alignment_base_import_does_not_import_numpy_or_biopython_align(self):
        code = """
import sys
import phykit.services.alignment.base as module
assert hasattr(module.np, "__getattr__")
assert callable(module.get_alignment_and_format_helper)
assert "typing" not in sys.modules
assert "numpy" not in sys.modules
assert "Bio.Align" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
assert "phykit.helpers.files" not in sys.modules
"""
        subprocess.run([sys.executable, "-c", code], check=True)

    def test_alignment_base_lazy_numpy_caches_resolved_attributes(self):
        lazy_np = alignment_base_module._LazyNumpy()

        zeros_attr = lazy_np.zeros

        assert lazy_np.__dict__["zeros"] is zeros_attr
        assert lazy_np.zeros is zeros_attr
        assert lazy_np._module is not None

    def test_alignment_json_output_modules_do_not_import_heavy_helpers(self):
        code = """
import importlib
import sys

modules = [
    "phykit.services.alignment.taxon_groups",
    "phykit.services.alignment.dna_threader",
    "phykit.services.alignment.occupancy_filter",
    "phykit.services.alignment.rcv",
]

for module_name in modules:
    module = importlib.import_module(module_name)
    assert callable(module.print_json)

assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "Bio.SeqIO" not in sys.modules
assert "Bio.SeqIO.FastaIO" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "numpy" not in sys.modules
"""
        subprocess.run([sys.executable, "-c", code], check=True)

    def test_init_sets_alignment_file_path(self, args):
        rcv = RelativeCompositionVariability(args)
        assert rcv.alignment_file_path == args.alignment
        assert rcv.output_file_path is None

    def test_relative_composition_variability(self, mocker, alignment_simple, args):
        mocker.patch(
            "phykit.services.alignment.rcv.RelativeCompositionVariability.get_alignment_and_format",
            return_value=(alignment_simple, "fa", True),
        )
        rcv = RelativeCompositionVariability(args)
        relative_composition_variability = rcv.calculate_rcv()
        assert isinstance(relative_composition_variability, float)
        assert isclose(relative_composition_variability, 0.292, rel_tol=0.001)

    def test_single_record_rcv_returns_before_sequence_materialization(
        self, mocker, args
    ):
        class UnstringableSequence:
            def __str__(self):
                raise AssertionError("single-record RCV should not inspect sequence")

        alignment = [SimpleNamespace(seq=UnstringableSequence(), id="t1")]
        rcv = RelativeCompositionVariability(args)
        mocker.patch.object(
            RelativeCompositionVariability,
            "get_alignment_and_format",
            return_value=(alignment, "fasta", False),
        )

        assert rcv.calculate_rcv() == 0.0

    def test_relative_composition_variability_handles_lowercase_ambiguity(self, mocker, args):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AnXT"), id="t1"),
                SeqRecord(Seq("ACxT"), id="t2"),
                SeqRecord(Seq("AG-T"), id="t3"),
            ]
        )
        rcv = RelativeCompositionVariability(args)

        mocked_get_alignment = mocker.patch.object(
            RelativeCompositionVariability,
            "get_alignment_and_format",
            return_value=(alignment, "fasta", False),
        )
        assert isclose(rcv.calculate_rcv(), 0.3333, rel_tol=0.001)

        mocked_get_alignment.return_value = (alignment, "fasta", True)
        assert isclose(rcv.calculate_rcv(), 0.4444, rel_tol=0.001)

    def test_relative_composition_variability_ascii_path_uses_lookup(self, mocker, args):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AnXT"), id="t1"),
                SeqRecord(Seq("ACxT"), id="t2"),
                SeqRecord(Seq("AG-T"), id="t3"),
            ]
        )
        rcv = RelativeCompositionVariability(args)
        mocker.patch.object(
            RelativeCompositionVariability,
            "get_alignment_and_format",
            return_value=(alignment, "fasta", False),
        )
        mocker.patch(
            "phykit.services.alignment.base.np.isin",
            side_effect=AssertionError("ASCII path should use the lookup table"),
        )

        assert isclose(rcv.calculate_rcv(), 0.3333, rel_tol=0.001)

    def test_relative_composition_variability_counts_valid_lengths_with_count_nonzero(
        self, mocker, args
    ):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AnXT"), id="t1"),
                SeqRecord(Seq("ACxT"), id="t2"),
                SeqRecord(Seq("AG-T"), id="t3"),
            ]
        )
        rcv = RelativeCompositionVariability(args)
        mocker.patch.object(
            RelativeCompositionVariability,
            "get_alignment_and_format",
            return_value=(alignment, "fasta", False),
        )
        count_nonzero_spy = mocker.spy(alignment_base_module.np, "count_nonzero")

        assert isclose(rcv.calculate_rcv(), 0.3333, rel_tol=0.001)
        assert any(
            call.kwargs.get("axis") == 1
            for call in count_nonzero_spy.call_args_list
        )

    def test_relative_composition_variability_final_total_uses_array_sum(
        self, mocker, args
    ):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGT"), id="t1"),
                SeqRecord(Seq("AGGT"), id="t2"),
            ]
        )
        rcv = RelativeCompositionVariability(args)
        mocker.patch.object(
            RelativeCompositionVariability,
            "get_alignment_and_format",
            return_value=(alignment, "fasta", True),
        )
        sum_spy = mocker.spy(alignment_base_module.np, "sum")

        assert rcv.calculate_rcv() > 0.0
        assert all("axis" in call.kwargs for call in sum_spy.call_args_list)

    def test_rcv_row_sums_use_array_reduction_for_narrow_matrices(self, monkeypatch):
        abs_diffs = alignment_base_module.np.array(
            [
                [0.0, 1.0, 0.5, 0.5],
                [0.25, 0.25, 0.25, 0.25],
            ],
            dtype=alignment_base_module.np.float64,
        )
        expected = abs_diffs.sum(axis=1)

        def fail_sum(*_args, **_kwargs):
            raise AssertionError("narrow RCV row sums should use ndarray.sum")

        monkeypatch.setattr(alignment_base_module.np, "sum", fail_sum)

        observed = alignment_base_module._rcv_row_sums(abs_diffs)

        alignment_base_module.np.testing.assert_allclose(observed, expected)

    def test_rcv_column_totals_use_array_reduction_for_common_matrices(
        self, monkeypatch
    ):
        count_matrix = alignment_base_module.np.array(
            [
                [10.0, 2.0, 5.0, 1.0],
                [4.0, 8.0, 3.0, 7.0],
                [6.0, 1.0, 9.0, 2.0],
            ],
            dtype=alignment_base_module.np.float64,
        )
        expected = count_matrix.sum(axis=0)

        def fail_sum(*_args, **_kwargs):
            raise AssertionError("common RCV column totals should use ndarray.sum")

        monkeypatch.setattr(alignment_base_module.np, "sum", fail_sum)

        observed = alignment_base_module._rcv_column_totals(count_matrix)

        alignment_base_module.np.testing.assert_allclose(observed, expected)

    def test_rcv_column_totals_preserve_large_matrix_sum_path(self, monkeypatch):
        count_matrix = alignment_base_module.np.ones(
            (1201, 20),
            dtype=alignment_base_module.np.float64,
        )
        original_sum = alignment_base_module.np.sum
        calls = []

        def sum_spy(values, *args, **kwargs):
            calls.append((values.shape, kwargs.get("axis")))
            return original_sum(values, *args, **kwargs)

        monkeypatch.setattr(alignment_base_module.np, "sum", sum_spy)

        observed = alignment_base_module._rcv_column_totals(count_matrix)

        alignment_base_module.np.testing.assert_allclose(
            observed,
            alignment_base_module.np.full(20, 1201.0),
        )
        assert calls == [((1201, 20), 0)]

    def test_relative_composition_variability_protein_ascii_uses_bincount(
        self, mocker, args
    ):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACDEFGHIKLMNPQRSTVWY"), id="t1"),
                SeqRecord(Seq("ACDEFGHIKLMNPQRSTVWA"), id="t2"),
            ]
        )
        rcv = RelativeCompositionVariability(args)
        mocker.patch.object(
            RelativeCompositionVariability,
            "get_alignment_and_format",
            return_value=(alignment, "fasta", True),
        )
        bincount_spy = mocker.spy(alignment_base_module.np, "bincount")

        assert rcv.calculate_rcv() > 0.0
        assert bincount_spy.call_count == len(alignment)

    def test_ascii_rcv_count_matrix_large_short_clean_uses_single_bincount(
        self,
        mocker,
    ):
        alphabet = b"ACDEFGHIKLMNPQRSTVWY"
        matrix = alignment_base_module.np.tile(
            alignment_base_module.np.frombuffer(
                alphabet,
                dtype=alignment_base_module.np.uint8,
            ),
            (10_000, 1),
        )
        unique_chars = alignment_base_module.np.unique(matrix)
        bincount_spy = mocker.spy(alignment_base_module.np, "bincount")

        observed = alignment_base_module._ascii_rcv_count_matrix(
            matrix,
            unique_chars,
            None,
        )

        expected = alignment_base_module.np.ones(
            (10_000, len(unique_chars)),
            dtype=alignment_base_module.np.float64,
        )
        alignment_base_module.np.testing.assert_array_equal(observed, expected)
        assert bincount_spy.call_count == 1

    def test_ascii_rcv_count_matrix_default_threshold_uses_single_bincount(
        self,
        mocker,
    ):
        alphabet = b"ACDEFGHIKLMNPQRSTVWY"
        matrix = alignment_base_module.np.resize(
            alignment_base_module.np.frombuffer(
                alphabet,
                dtype=alignment_base_module.np.uint8,
            ),
            2_048 * 512,
        ).reshape(2_048, 512)
        unique_chars = alignment_base_module.np.unique(matrix)
        bincount_spy = mocker.spy(alignment_base_module.np, "bincount")

        observed = alignment_base_module._ascii_rcv_count_matrix(
            matrix,
            unique_chars,
            None,
        )

        expected = alignment_base_module.np.zeros(
            (matrix.shape[0], len(unique_chars)),
            dtype=alignment_base_module.np.float64,
        )
        for row_idx, row in enumerate(matrix):
            expected[row_idx] = alignment_base_module.np.bincount(
                row,
                minlength=256,
            )[unique_chars]
        alignment_base_module.np.testing.assert_array_equal(observed, expected)
        assert bincount_spy.call_count == matrix.shape[0] + 1

    def test_ascii_rcv_count_matrix_large_short_gappy_uses_single_bincount(
        self,
        mocker,
    ):
        alphabet = b"ACGTN-?X*"
        matrix = alignment_base_module.np.tile(
            alignment_base_module.np.frombuffer(
                alphabet,
                dtype=alignment_base_module.np.uint8,
            ),
            (10_000, 1),
        )
        invalid_lookup = alignment_base_module._get_invalid_lookup(is_protein=False)
        observed_chars = alignment_base_module.np.unique(matrix)
        unique_chars = observed_chars[~invalid_lookup[observed_chars]]
        valid_mask = ~invalid_lookup[matrix]
        bincount_spy = mocker.spy(alignment_base_module.np, "bincount")

        observed = alignment_base_module._ascii_rcv_count_matrix(
            matrix,
            unique_chars,
            valid_mask,
        )

        expected = alignment_base_module.np.ones(
            (10_000, len(unique_chars)),
            dtype=alignment_base_module.np.float64,
        )
        alignment_base_module.np.testing.assert_array_equal(observed, expected)
        assert bincount_spy.call_count == 1

    def test_relative_composition_variability_no_gap_ascii_uses_full_lengths(
        self, mocker, args
    ):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGT"), id="t1"),
                SeqRecord(Seq("AGGT"), id="t2"),
            ]
        )
        rcv = RelativeCompositionVariability(args)
        mocker.patch.object(
            RelativeCompositionVariability,
            "get_alignment_and_format",
            return_value=(alignment, "fasta", True),
        )
        full_spy = mocker.spy(alignment_base_module.np, "full")
        mocker.patch(
            "phykit.services.alignment.base.np.isin",
            side_effect=AssertionError("ASCII path should not use Unicode isin"),
        )

        assert rcv.calculate_rcv() > 0.0
        full_spy.assert_called_once_with(
            len(alignment),
            alignment.get_alignment_length(),
            dtype=alignment_base_module.np.float64,
        )

    def test_clean_nucleotide_rcv_skips_full_symbol_discovery(
        self,
        monkeypatch,
        args,
    ):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGT"), id="t1"),
                SeqRecord(Seq("AGGT"), id="t2"),
                SeqRecord(Seq("TCGA"), id="t3"),
            ]
        )
        rcv = RelativeCompositionVariability(args)
        monkeypatch.setattr(
            rcv,
            "get_alignment_and_format",
            lambda: (alignment, "fasta", False),
        )

        def fail_unique(*args, **kwargs):
            raise AssertionError("clean nucleotides should use known symbol codes")

        monkeypatch.setattr(alignment_base_module.np, "unique", fail_unique)

        assert rcv.calculate_rcv() == pytest.approx(2 / 9)

    def test_relative_composition_variability_identical_sequences_skip_matrix(
        self, mocker, args
    ):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGTACGT"), id="t1"),
                SeqRecord(Seq("acgtacgt"), id="t2"),
                SeqRecord(Seq("ACGTACGT"), id="t3"),
            ]
        )
        rcv = RelativeCompositionVariability(args)
        mocker.patch.object(
            RelativeCompositionVariability,
            "get_alignment_and_format",
            return_value=(alignment, "fasta", False),
        )
        mocker.patch(
            "phykit.services.alignment.base.np.frombuffer",
            side_effect=AssertionError("identical sequences should skip matrix work"),
        )

        assert rcv.calculate_rcv() == 0.0

    def test_identical_sequence_helper_does_not_slice_rows(self):
        class NoSliceList(list):
            def __getitem__(self, key):
                if isinstance(key, slice):
                    raise AssertionError("identical-sequence scan should not slice")
                return super().__getitem__(key)

        sequences = NoSliceList(["ACGT", "ACGT", "ACGT"])

        assert alignment_base_module._all_sequences_identical(sequences) is True

    def test_relative_composition_variability_unicode_fallback(self, mocker, args):
        class DummyAlignment(list):
            def get_alignment_length(self):
                return 2

        alignment = DummyAlignment(
            [
                SimpleNamespace(seq="A\u00d1", id="t1"),
                SimpleNamespace(seq="A\u00d1", id="t2"),
                SimpleNamespace(seq="T\u00d1", id="t3"),
            ]
        )
        rcv = RelativeCompositionVariability(args)
        mocker.patch.object(
            RelativeCompositionVariability,
            "get_alignment_and_format",
            return_value=(alignment, "fasta", False),
        )

        assert isclose(rcv.calculate_rcv(), 0.4444, rel_tol=0.001)

    def test_process_args_defaults_json_false(self):
        parsed = RelativeCompositionVariability(Namespace(alignment="x.fa")).process_args(
            Namespace(alignment="x.fa")
        )
        assert parsed["json_output"] is False

    def test_run_prints_rcv(self, mocker, capsys):
        rcv = RelativeCompositionVariability(Namespace(alignment="x.fa", json=False))
        mocker.patch.object(RelativeCompositionVariability, "calculate_rcv", return_value=0.123456)
        rcv.run()
        out, _ = capsys.readouterr()
        assert out.strip() == "0.1235"

    def test_run_json_output(self, mocker):
        rcv = RelativeCompositionVariability(Namespace(alignment="x.fa", json=True))
        mocker.patch.object(RelativeCompositionVariability, "calculate_rcv", return_value=0.56789)
        mocked_json = mocker.patch.object(rcv_module, "print_json")
        rcv.run()
        mocked_json.assert_called_once_with({"rcv": 0.5679})
