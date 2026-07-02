import pytest
import subprocess
import sys
from argparse import Namespace
from types import SimpleNamespace

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.mask_alignment import MaskAlignment
import phykit.services.alignment.mask_alignment as mask_alignment_module


@pytest.fixture
def args():
    kwargs = dict(
        alignment="/some/path/to/file.fa",
        max_gap=1.0,
        min_occupancy=0.0,
        max_entropy=None,
    )
    return Namespace(**kwargs)


class TestMaskAlignment(object):
    def test_module_import_does_not_import_numpy_or_biopython_align(self):
        code = """
import sys
import phykit.services.alignment.mask_alignment as module
assert hasattr(module.np, "__getattr__")
assert "typing" not in sys.modules
assert "numpy" not in sys.modules
assert "Bio.Align" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
"""
        subprocess.run([sys.executable, "-c", code], check=True)

    def test_init_sets_alignment_file_path(self, args):
        masker = MaskAlignment(args)
        assert masker.alignment_file_path == args.alignment

    def test_keep_mask_by_gap_and_occupancy(self, alignment_simple, args):
        args.max_gap = 0.3
        args.min_occupancy = 0.8
        masker = MaskAlignment(args)
        keep_mask = masker.calculate_keep_mask(alignment_simple, is_protein=False)
        assert keep_mask.tolist() == [True, False, True, False, True, True]

    def test_keep_mask_by_entropy(self, alignment_simple, args):
        args.max_entropy = 0.5
        masker = MaskAlignment(args)
        keep_mask = masker.calculate_keep_mask(alignment_simple, is_protein=False)
        assert keep_mask.tolist() == [True, False, False, True, False, False]

    def test_keep_mask_ascii_entropy_counts_with_count_nonzero(self, mocker, args):
        args.max_entropy = 1.0
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGT"), id="a"),
                SeqRecord(Seq("ATGT"), id="b"),
                SeqRecord(Seq("TCGT"), id="c"),
            ]
        )
        masker = MaskAlignment(args)
        count_nonzero_spy = mocker.spy(mask_alignment_module.np, "count_nonzero")

        keep_mask = masker.calculate_keep_mask(alignment, is_protein=False)

        assert keep_mask.tolist() == [True, True, True, True]
        assert any(
            call.kwargs.get("axis") == 0
            for call in count_nonzero_spy.call_args_list
        )

    def test_keep_mask_by_entropy_handles_lowercase_and_all_invalid(self, args):
        args.max_entropy = 0.0
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("Ac-"), id="a"),
                SeqRecord(Seq("AgN"), id="b"),
                SeqRecord(Seq("Acx"), id="c"),
            ]
        )
        masker = MaskAlignment(args)

        keep_mask = masker.calculate_keep_mask(alignment, is_protein=False)

        assert keep_mask.tolist() == [True, False, True]

    def test_keep_mask_ascii_path_uses_lookup(self, mocker, args):
        args.max_entropy = 0.0
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("Ac-"), id="a"),
                SeqRecord(Seq("AgN"), id="b"),
                SeqRecord(Seq("Acx"), id="c"),
            ]
        )
        masker = MaskAlignment(args)
        mocker.patch(
            "phykit.services.alignment.mask_alignment.np.isin",
            side_effect=AssertionError("ASCII path should use the lookup table"),
        )

        keep_mask = masker.calculate_keep_mask(alignment, is_protein=False)

        assert keep_mask.tolist() == [True, False, True]

    def test_keep_mask_all_pass_thresholds_short_circuit(self, mocker, args):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("A-GN"), id="a"),
                SeqRecord(Seq("T?CX"), id="b"),
            ]
        )
        masker = MaskAlignment(args)
        mocker.patch(
            "phykit.services.alignment.mask_alignment.np.mean",
            side_effect=AssertionError("all-pass thresholds should skip occupancy mean"),
        )

        keep_mask = masker.calculate_keep_mask(alignment, is_protein=False)

        assert keep_mask.tolist() == [True, True, True, True]

    def test_keep_mask_all_pass_thresholds_skip_sequence_materialization(
        self, args
    ):
        class UnstringableSequence:
            def __str__(self):
                raise AssertionError("all-pass mask should not inspect sequences")

        class DummyAlignment(list):
            def get_alignment_length(self):
                return 4

        alignment = DummyAlignment([SimpleNamespace(seq=UnstringableSequence())])
        masker = MaskAlignment(args)

        keep_mask = masker.calculate_keep_mask(alignment, is_protein=False)

        assert keep_mask.tolist() == [True, True, True, True]

    def test_keep_mask_clean_entropy_skips_occupancy(self, mocker, args):
        args.max_entropy = 0.0
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AAAA"), id="a"),
                SeqRecord(Seq("AAAA"), id="b"),
            ]
        )
        masker = MaskAlignment(args)
        mean_spy = mocker.spy(mask_alignment_module.np, "mean")
        mocker.patch(
            "phykit.services.alignment.mask_alignment.np.sum",
            side_effect=AssertionError(
                "single-symbol entropy path should skip counts"
            ),
        )

        keep_mask = masker.calculate_keep_mask(alignment, is_protein=False)

        assert keep_mask.tolist() == [True, True, True, True]
        assert mean_spy.call_count == 0

    def test_keep_mask_identical_sequences_uses_byte_mask(self, mocker, args):
        args.max_gap = 0.5
        args.min_occupancy = 0.5
        args.max_entropy = 0.0
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AcGN-?"), id="a"),
                SeqRecord(Seq("aCgN-?"), id="b"),
                SeqRecord(Seq("ACGN-?"), id="c"),
            ]
        )
        masker = MaskAlignment(args)
        frombuffer_spy = mocker.spy(mask_alignment_module.np, "frombuffer")
        fromiter_spy = mocker.spy(mask_alignment_module.np, "fromiter")
        mocker.patch(
            "phykit.services.alignment.mask_alignment.np.mean",
            side_effect=AssertionError("identical alignments should skip occupancy mean"),
        )

        keep_mask = masker.calculate_keep_mask(alignment, is_protein=False)

        assert keep_mask.tolist() == [True, True, True, False, False, False]
        assert frombuffer_spy.call_count == 1
        assert fromiter_spy.call_count == 0

    def test_identical_sequence_helper_does_not_slice_rows(self):
        class NoSliceList(list):
            def __getitem__(self, key):
                if isinstance(key, slice):
                    raise AssertionError("identical-sequence scan should not slice")
                return super().__getitem__(key)

        sequences = NoSliceList(["ACGT", "ACGT", "ACGT"])

        assert mask_alignment_module._all_sequences_identical(sequences) is True

    def test_keep_mask_clean_no_entropy_skips_matrix_and_occupancy(self, mocker, args):
        args.max_gap = 0.2
        args.min_occupancy = 0.8
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("acgt"), id="a"),
                SeqRecord(Seq("tgca"), id="b"),
            ]
        )
        masker = MaskAlignment(args)
        mocker.patch(
            "phykit.services.alignment.mask_alignment.np.frombuffer",
            side_effect=AssertionError("clean no-entropy path should skip matrix setup"),
        )
        mocker.patch(
            "phykit.services.alignment.mask_alignment.np.mean",
            side_effect=AssertionError("clean no-entropy path should skip occupancy"),
        )

        keep_mask = masker.calculate_keep_mask(alignment, is_protein=False)

        assert keep_mask.tolist() == [True, True, True, True]

    def test_keep_mask_protein_entropy_uses_block_counts(self, mocker, args):
        args.max_entropy = 0.1
        alphabet = "ACDEFGHIKL"
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq(symbol + "A"), id=f"t{idx}")
                for idx, symbol in enumerate(alphabet)
            ]
        )
        masker = MaskAlignment(args)
        mocker.patch(
            "phykit.services.alignment.mask_alignment.np.isin",
            side_effect=AssertionError("ASCII protein path should use byte counts"),
        )
        mocker.patch(
            "phykit.services.alignment.mask_alignment.np.mean",
            side_effect=AssertionError("clean entropy path should skip occupancy"),
        )
        entropy_blocks = mocker.spy(
            mask_alignment_module,
            "_column_entropies_from_ascii_codes",
        )

        keep_mask = masker.calculate_keep_mask(alignment, is_protein=True)

        assert entropy_blocks.call_args.args[1] is None
        assert keep_mask.tolist() == [False, True]

    def test_column_entropies_from_ascii_codes_avoids_tiled_offsets(
        self, mocker
    ):
        alignment_array = mask_alignment_module.np.array(
            [
                list(b"ACDEFX"),
                list(b"ACGHIX"),
                list(b"LMNPQ-"),
            ],
            dtype=mask_alignment_module.np.uint8,
        )
        invalid_lookup = mask_alignment_module._get_invalid_lookup(is_protein=True)
        valid_mask = ~invalid_lookup[alignment_array]
        valid_symbols = mask_alignment_module.np.unique(alignment_array[valid_mask])
        mocker.patch(
            "phykit.services.alignment.mask_alignment.np.tile",
            side_effect=AssertionError("ASCII entropy path should avoid tiled offsets"),
        )

        entropies = mask_alignment_module._column_entropies_from_ascii_codes(
            alignment_array,
            valid_mask,
            valid_symbols,
            block_size=2,
        )

        assert entropies.shape == (alignment_array.shape[1],)
        assert entropies[-1] == 0.0

    def test_column_entropies_from_ascii_codes_uses_masked_log2(
        self, monkeypatch
    ):
        alignment_array = mask_alignment_module.np.array(
            [
                list(b"ACDE"),
                list(b"AAAA"),
                list(b"AAAA"),
            ],
            dtype=mask_alignment_module.np.uint8,
        )
        valid_symbols = mask_alignment_module.np.unique(alignment_array)
        real_log2 = mask_alignment_module.np.log2
        log_calls = []

        def checked_log2(values, *args, **kwargs):
            assert "out" in kwargs
            assert "where" in kwargs
            log_calls.append(values.shape)
            return real_log2(values, *args, **kwargs)

        monkeypatch.setattr(
            mask_alignment_module.np,
            "log2",
            checked_log2,
            raising=False,
        )

        entropies = mask_alignment_module._column_entropies_from_ascii_codes(
            alignment_array,
            None,
            valid_symbols,
            block_size=2,
        )

        assert log_calls
        assert entropies.shape == (alignment_array.shape[1],)

    def test_entropy_columns_from_probabilities_matches_explicit_sum(self):
        np = mask_alignment_module.np
        probs = np.array(
            [
                [0.2, 0.0, 0.5],
                [0.3, 1.0, 0.5],
                [0.5, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
            ],
            dtype=np.float64,
        )
        log_probs = np.zeros_like(probs)
        positive = probs > 0
        np.log2(probs, out=log_probs, where=positive)

        observed = mask_alignment_module._entropy_columns_from_probabilities(
            probs.copy(),
            log_probs,
        )
        expected = -np.sum(probs * log_probs, axis=0)

        np.testing.assert_allclose(observed, expected)

    def test_entropy_columns_from_probabilities_uses_column_dot(self, monkeypatch):
        np = mask_alignment_module.np
        probs = np.array(
            [
                [0.2, 0.0, 0.5],
                [0.3, 1.0, 0.5],
                [0.5, 0.0, 0.0],
            ],
            dtype=np.float64,
        )
        log_probs = np.zeros_like(probs)
        positive = probs > 0
        np.log2(probs, out=log_probs, where=positive)
        expected = -np.einsum("ij,ij->j", probs, log_probs)

        def fail_sum(*_args, **_kwargs):
            raise AssertionError("entropy columns should use a column dot")

        monkeypatch.setattr(mask_alignment_module.np, "sum", fail_sum)

        observed = mask_alignment_module._entropy_columns_from_probabilities(
            probs.copy(),
            log_probs,
        )

        np.testing.assert_allclose(observed, expected)

    def test_keep_mask_unicode_fallback(self, args):
        class DummyAlignment(list):
            def get_alignment_length(self):
                return 2

        args.max_entropy = 1.0
        alignment = DummyAlignment(
            [
                SimpleNamespace(seq="A\u00d1", id="a"),
                SimpleNamespace(seq="A\u00d1", id="b"),
                SimpleNamespace(seq="T\u00d1", id="c"),
            ]
        )
        masker = MaskAlignment(args)

        keep_mask = masker.calculate_keep_mask(alignment, is_protein=False)

        assert keep_mask.tolist() == [True, True]

    def test_process_args_json_default(self):
        parsed = MaskAlignment(
            Namespace(alignment="x.fa", max_gap=0.8, min_occupancy=0.2, max_entropy=None)
        ).process_args(
            Namespace(alignment="x.fa", max_gap=0.8, min_occupancy=0.2, max_entropy=None)
        )
        assert parsed["json_output"] is False

    def test_validate_thresholds_invalid_values_exit(self, capsys):
        with pytest.raises(SystemExit) as exc:
            MaskAlignment(
                Namespace(alignment="x.fa", max_gap=1.2, min_occupancy=0.0, max_entropy=None)
            )._validate_thresholds()
        assert exc.value.code == 2
        out, _ = capsys.readouterr()
        assert "max_gap must be between 0 and 1." in out

        with pytest.raises(SystemExit) as exc:
            MaskAlignment(
                Namespace(alignment="x.fa", max_gap=0.8, min_occupancy=-0.1, max_entropy=None)
            )._validate_thresholds()
        assert exc.value.code == 2
        out, _ = capsys.readouterr()
        assert "min_occupancy must be between 0 and 1." in out

        with pytest.raises(SystemExit) as exc:
            MaskAlignment(
                Namespace(alignment="x.fa", max_gap=0.8, min_occupancy=0.1, max_entropy=-0.5)
            )._validate_thresholds()
        assert exc.value.code == 2
        out, _ = capsys.readouterr()
        assert "max_entropy must be >= 0." in out

    def test_apply_mask(self, alignment_simple, args):
        masker = MaskAlignment(args)
        masked = masker.apply_mask(alignment_simple, keep_mask=mask_alignment_module.np.array([True, False, True, False, True, False]))
        assert set(masked.keys()) == {record.id for record in alignment_simple}
        assert all(len(seq) == 3 for seq in masked.values())

    def test_apply_mask_uppercases_ascii_sequences(self, args):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("acgt"), id="a"),
                SeqRecord(Seq("tgca"), id="b"),
            ]
        )
        masker = MaskAlignment(args)

        masked = masker.apply_mask(
            alignment,
            keep_mask=mask_alignment_module.np.array([True, False, True, False]),
        )

        assert masked == {"a": "AG", "b": "TC"}

    def test_apply_mask_all_keep_returns_uppercase_sequences_directly(
        self, mocker, args
    ):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("acgt"), id="a"),
                SeqRecord(Seq("tgca"), id="b"),
            ]
        )
        masker = MaskAlignment(args)
        mocker.patch(
            "phykit.services.alignment.mask_alignment.np.frombuffer",
            side_effect=AssertionError("all-keep mask should not build a matrix"),
        )

        masked = masker.apply_mask(
            alignment,
            keep_mask=mask_alignment_module.np.array([True, True, True, True]),
        )

        assert masked == {"a": "ACGT", "b": "TGCA"}

    def test_apply_mask_all_trim_returns_empty_sequences_directly(
        self, mocker, args
    ):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("acgt"), id="a"),
                SeqRecord(Seq("tgca"), id="b"),
            ]
        )
        masker = MaskAlignment(args)
        mocker.patch(
            "phykit.services.alignment.mask_alignment.np.frombuffer",
            side_effect=AssertionError("all-trim mask should not build a matrix"),
        )

        masked = masker.apply_mask(
            alignment,
            keep_mask=mask_alignment_module.np.array([False, False, False, False]),
        )

        assert masked == {"a": "", "b": ""}

    def test_apply_mask_all_trim_mismatched_length_uses_matrix_path(
        self, mocker, args
    ):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("acgt"), id="a"),
                SeqRecord(Seq("tgca"), id="b"),
            ]
        )
        masker = MaskAlignment(args)
        frombuffer_spy = mocker.spy(mask_alignment_module.np, "frombuffer")

        with pytest.raises(ValueError):
            masker.apply_mask(
                alignment,
                keep_mask=mask_alignment_module.np.array([False, False]),
            )

        assert frombuffer_spy.call_count == 1

    def test_apply_mask_falls_back_for_non_ascii_sequences(self, args):
        alignment = [
            SimpleNamespace(seq="A\u00d1GT", id="a"),
            SimpleNamespace(seq="TCGA", id="b"),
        ]
        masker = MaskAlignment(args)

        masked = masker.apply_mask(
            alignment,
            keep_mask=mask_alignment_module.np.array([True, True, False, False]),
        )

        assert masked == {"a": "A\u00d1", "b": "TC"}

    def test_run_json_output(self, alignment_simple, mocker):
        masker = MaskAlignment(
            Namespace(alignment="x.fa", max_gap=1.0, min_occupancy=0.0, max_entropy=None, json=True)
        )
        mocker.patch.object(MaskAlignment, "get_alignment_and_format", return_value=(alignment_simple, "fasta", False))
        mocker.patch.object(MaskAlignment, "calculate_keep_mask", return_value=mask_alignment_module.np.array([True, True, False, False, True, True]))
        mocker.patch.object(MaskAlignment, "apply_mask", return_value={"a": "ACGT", "b": "ACGT"})
        mocked_json = mocker.patch("phykit.services.alignment.mask_alignment.print_json")

        masker.run()

        payload = mocked_json.call_args.args[0]
        assert payload["kept_sites"] == 4
        assert payload["total_sites"] == 6
        assert payload["rows"] == payload["taxa"]

    def test_run_prints_fasta(self, alignment_simple, capsys, mocker):
        masker = MaskAlignment(
            Namespace(alignment="x.fa", max_gap=1.0, min_occupancy=0.0, max_entropy=None, json=False)
        )
        mocker.patch.object(MaskAlignment, "get_alignment_and_format", return_value=(alignment_simple, "fasta", False))
        mocker.patch.object(MaskAlignment, "calculate_keep_mask", return_value=mask_alignment_module.np.array([True, False, True, False, True, False]))
        mocker.patch.object(MaskAlignment, "apply_mask", return_value={"a": "AAA", "b": "CCC"})

        masker.run()

        out, _ = capsys.readouterr()
        assert out == ">a\nAAA\n>b\nCCC\n"
