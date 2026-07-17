from argparse import Namespace
import subprocess
import sys
from types import SimpleNamespace

import pytest
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.gc_content import GCContent
import phykit.services.alignment.gc_content as gc_content_module


def _alignment(records):
    return MultipleSeqAlignment(records)


@pytest.fixture
def args():
    return Namespace(fasta="/some/path/to/file.fa", verbose=False)


class TestGCContent:
    def test_module_import_does_not_import_numpy_or_biopython_align(self):
        code = """
import sys
import phykit.services.alignment.gc_content as module
assert hasattr(module.np, "__getattr__")
assert callable(module.get_alignment_and_format)
assert "typing" not in sys.modules
assert "numpy" not in sys.modules
assert "Bio.Align" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
assert "phykit.helpers.files" not in sys.modules
"""
        subprocess.run([sys.executable, "-c", code], check=True)

    def test_lazy_numpy_proxy_caches_resolved_attributes(self):
        lazy_np = gc_content_module._LazyNumpy()

        first = lazy_np.count_nonzero
        second = lazy_np.count_nonzero

        assert first is second
        assert lazy_np.__dict__["count_nonzero"] is first

    def test_init_sets_expected_attrs(self, args):
        service = GCContent(args)
        assert service.fasta == args.fasta
        assert service.verbose is False
        assert service.json_output is False

    def test_common_upper_sequence_matches_case_insensitive_rows(self):
        assert (
            gc_content_module._common_upper_sequence(["acgt", "ACGT", "acgt"])
            == "ACGT"
        )
        assert gc_content_module._common_upper_sequence(["acgt", "ACGA"]) is None

    def test_gc_counts_from_upper_sequence_counts_ascii_bytes_and_unicode(self):
        assert gc_content_module._gc_counts_from_upper_sequence(
            "ACGTGCNN--??XX",
            is_protein=False,
        ) == (6, 4)
        assert gc_content_module._gc_counts_from_upper_sequence(
            "ACGTGCNN--??XX",
            is_protein=True,
        ) == (8, 4)
        assert gc_content_module._gc_counts_from_upper_sequence(
            "G\u03a9N-",
            is_protein=False,
        ) == (2, 1)

    def test_gc_counts_from_mixed_sequence_avoids_upper_copy(self):
        class NoUpperString(str):
            def upper(self):
                raise AssertionError(
                    "Unicode GC fallback should count mixed-case symbols directly"
                )

        sequence = NoUpperString("gC\u03a9nX-")

        assert gc_content_module._gc_counts_from_mixed_sequence(
            sequence,
            is_protein=False,
        ) == (3, 2)
        assert gc_content_module._gc_counts_from_mixed_sequence(
            sequence,
            is_protein=True,
        ) == (4, 2)

    def test_invalid_byte_scan_short_buffers_avoid_slice_checks(self):
        class NoSliceBytes(bytes):
            def __getitem__(self, key):
                if isinstance(key, slice):
                    raise AssertionError("short invalid scan should not slice")
                return super().__getitem__(key)

        clean = NoSliceBytes(b"ACGT" * 10)
        gapped = NoSliceBytes(b"ACGTN")

        assert gc_content_module._has_invalid_bytes(
            clean,
            gc_content_module._DNA_INVALID_BYTES,
        ) is False
        assert gc_content_module._has_invalid_bytes(
            gapped,
            gc_content_module._DNA_INVALID_BYTES,
        ) is True

    def test_calculate_gc_total_value(self, args):
        service = GCContent(args)
        records = _alignment(
            [
                SeqRecord(Seq("ACGT"), id="a"),
                SeqRecord(Seq("GGTT"), id="b"),
            ]
        )
        assert service.calculate_gc_total_value(records, is_protein=False) == 0.5

    def test_calculate_gc_total_value_ignores_ambiguous_sites_and_uppercases(self, args):
        service = GCContent(args)
        records = _alignment(
            [
                SeqRecord(Seq("acgn"), id="a"),
                SeqRecord(Seq("G-C?"), id="b"),
            ]
        )
        assert service.calculate_gc_total_value(records, is_protein=False) == round(4 / 5, 4)

    def test_calculate_gc_total_value_exits_for_empty_cleaned_seq(self, mocker, args):
        service = GCContent(args)
        records = _alignment([SeqRecord(Seq("----"), id="a")])
        mocked_print = mocker.patch("builtins.print")
        with pytest.raises(SystemExit) as excinfo:
            service.calculate_gc_total_value(records, is_protein=False)
        assert excinfo.value.code == 2
        mocked_print.assert_called_once()

    def test_calculate_gc_per_sequence_data(self, args):
        service = GCContent(args)
        records = _alignment(
            [
                SeqRecord(Seq("G-C"), id="a"),  # cleaned GC = 2/2
                SeqRecord(Seq("A-T"), id="b"),  # cleaned GC = 0/2
            ]
        )
        assert service.calculate_gc_per_sequence_data(records, is_protein=False) == [
            ("a", 1.0),
            ("b", 0.0),
        ]

    def test_calculate_gc_per_sequence_data_ignores_ambiguous_sites_and_uppercases(self, args):
        service = GCContent(args)
        records = _alignment(
            [
                SeqRecord(Seq("gcn"), id="a"),
                SeqRecord(Seq("a-t"), id="b"),
            ]
        )
        assert service.calculate_gc_per_sequence_data(records, is_protein=False) == [
            ("a", 1.0),
            ("b", 0.0),
        ]

    def test_calculate_gc_per_sequence_data_uses_bounded_ascii_matrix_counts(
        self, args, mocker, monkeypatch
    ):
        service = GCContent(args)
        records = _alignment(
            [
                SeqRecord(Seq("gcN-"), id="a"),
                SeqRecord(Seq("AT??"), id="b"),
            ]
        )
        monkeypatch.setattr(gc_content_module, "_GC_PER_SEQUENCE_SCALAR_MAX_BYTES", 0)
        mocker.patch.object(
            gc_content_module.np,
            "count_nonzero",
            side_effect=AssertionError(
                "ASCII matrix path should use bounded row sums"
            ),
        )

        assert service.calculate_gc_per_sequence_data(records, is_protein=False) == [
            ("a", 1.0),
            ("b", 0.0),
        ]

    @pytest.mark.parametrize(
        ("sites", "expected_dtype"),
        [(0xFFFF, "uint16"), (0x10000, "uint32")],
    )
    def test_bounded_ascii_row_counts_preserve_site_count(
        self, sites, expected_dtype
    ):
        mask = gc_content_module.np.ones((1, sites), dtype=bool)

        counts = gc_content_module._bounded_ascii_row_counts(mask)

        assert int(counts[0]) == sites
        assert counts.dtype == gc_content_module.np.dtype(expected_dtype)

    def test_calculate_gc_per_sequence_data_no_gap_ascii_skips_valid_lookup(
        self, args, mocker, monkeypatch
    ):
        service = GCContent(args)
        records = _alignment(
            [
                SeqRecord(Seq("GCGT"), id="a"),
                SeqRecord(Seq("ATAT"), id="b"),
            ]
        )
        monkeypatch.setattr(gc_content_module, "_GC_PER_SEQUENCE_SCALAR_MAX_BYTES", 0)
        mocker.patch(
            "phykit.services.alignment.gc_content._get_valid_lookup",
            side_effect=AssertionError(
                "All-valid ASCII GC rows should not build valid lookup counts"
            ),
        )

        assert service.calculate_gc_per_sequence_data(records, is_protein=False) == [
            ("a", 0.75),
            ("b", 0.0),
        ]

    def test_calculate_gc_per_sequence_data_identical_sequences_skip_matrix(
        self, args, mocker, monkeypatch
    ):
        service = GCContent(args)
        records = _alignment(
            [
                SeqRecord(Seq("acgtGCNN--??XX"), id="a"),
                SeqRecord(Seq("ACGTGCNN--??XX"), id="b"),
                SeqRecord(Seq("acgtGCNN--??XX"), id="c"),
            ]
        )
        monkeypatch.setattr(gc_content_module, "_GC_PER_SEQUENCE_SCALAR_MAX_BYTES", 0)
        mocker.patch(
            "phykit.services.alignment.gc_content.np.frombuffer",
            side_effect=AssertionError(
                "identical GC rows should not build a byte matrix"
            ),
        )

        assert service.calculate_gc_per_sequence_data(records, is_protein=False) == [
            ("a", 4 / 6),
            ("b", 4 / 6),
            ("c", 4 / 6),
        ]

    def test_calculate_gc_per_sequence_data_variable_length_uses_byte_counts(
        self, args, mocker, monkeypatch
    ):
        service = GCContent(args)
        records = [
            SimpleNamespace(seq="GCN-", id="a"),
            SimpleNamespace(seq="A-C", id="b"),
        ]
        monkeypatch.setattr(gc_content_module, "_GC_PER_SEQUENCE_SCALAR_MAX_BYTES", 0)
        mocker.patch.object(
            gc_content_module.np,
            "frombuffer",
            side_effect=AssertionError(
                "variable-length ASCII fallback should use byte counts"
            ),
        )

        assert service.calculate_gc_per_sequence_data(records, is_protein=False) == [
            ("a", 1.0),
            ("b", 0.5),
        ]

    def test_calculate_gc_per_sequence_data_small_ascii_avoids_numpy(
        self, args, mocker
    ):
        service = GCContent(args)
        records = _alignment(
            [
                SeqRecord(Seq("gcN-"), id="a"),
                SeqRecord(Seq("AT??"), id="b"),
            ]
        )
        mocker.patch(
            "phykit.services.alignment.gc_content._gc_counts_from_ascii_matrix",
            side_effect=AssertionError(
                "small verbose GC rows should avoid matrix setup"
            ),
        )
        mocker.patch(
            "phykit.services.alignment.gc_content.np.frombuffer",
            side_effect=AssertionError(
                "small verbose GC rows should avoid NumPy setup"
            ),
        )

        assert service.calculate_gc_per_sequence_data(records, is_protein=False) == [
            ("a", 1.0),
            ("b", 0.0),
        ]

    def test_calculate_gc_total_value_uses_ascii_flat_path(self, args, mocker):
        service = GCContent(args)
        records = _alignment(
            [
                SeqRecord(Seq("gcN-"), id="a"),
                SeqRecord(Seq("AT??"), id="b"),
            ]
        )
        mocker.patch.object(
            gc_content_module.np,
            "sum",
            side_effect=AssertionError("ASCII total path should avoid row sums"),
        )

        assert service.calculate_gc_total_value(records, is_protein=False) == round(
            2 / 4,
            4,
        )

    def test_calculate_gc_total_value_small_ascii_avoids_numpy(
        self,
        args,
        mocker,
    ):
        service = GCContent(args)
        records = _alignment(
            [
                SeqRecord(Seq("GCGT"), id="a"),
                SeqRecord(Seq("ATAT"), id="b"),
            ]
        )
        mocker.patch(
            "phykit.services.alignment.gc_content.np.frombuffer",
            side_effect=AssertionError(
                "small ASCII total GC should use scalar byte counts"
            ),
        )

        assert service.calculate_gc_total_value(records, is_protein=False) == round(
            3 / 8,
            4,
        )

    def test_calculate_gc_total_value_invalid_ascii_uses_byte_counts(
        self, args, monkeypatch, mocker
    ):
        service = GCContent(args)
        records = _alignment(
            [
                SeqRecord(Seq("gcN-"), id="a"),
                SeqRecord(Seq("AT??"), id="b"),
            ]
        )
        monkeypatch.setattr(gc_content_module, "_GC_TOTAL_BYTE_COUNT_MIN_BYTES", 1)
        count_nonzero = mocker.patch(
            "phykit.services.alignment.gc_content.np.count_nonzero",
            side_effect=AssertionError(
                "invalid ASCII total GC should not use lookup reductions"
            ),
        )

        assert service.calculate_gc_total_value(records, is_protein=False) == round(
            2 / 4,
            4,
        )
        count_nonzero.assert_not_called()

    def test_calculate_gc_total_value_no_gap_ascii_skips_valid_lookup(
        self, args, mocker
    ):
        service = GCContent(args)
        records = _alignment(
            [
                SeqRecord(Seq("GCGT"), id="a"),
                SeqRecord(Seq("ATAT"), id="b"),
            ]
        )
        mocker.patch(
            "phykit.services.alignment.gc_content._get_valid_lookup",
            side_effect=AssertionError(
                "All-valid ASCII GC totals should not build valid lookup counts"
            ),
        )

        assert service.calculate_gc_total_value(records, is_protein=False) == round(
            3 / 8,
            4,
        )

    def test_calculate_gc_total_value_identical_sequences_skip_flat_matrix(
        self, args, mocker
    ):
        service = GCContent(args)
        records = _alignment(
            [
                SeqRecord(Seq("acgtGCNN--??XX"), id="a"),
                SeqRecord(Seq("ACGTGCNN--??XX"), id="b"),
                SeqRecord(Seq("acgtGCNN--??XX"), id="c"),
            ]
        )
        mocker.patch(
            "phykit.services.alignment.gc_content.np.frombuffer",
            side_effect=AssertionError(
                "identical GC totals should not build a flat byte matrix"
            ),
        )

        assert service.calculate_gc_total_value(records, is_protein=False) == round(
            4 / 6,
            4,
        )

    def test_calculate_gc_per_sequence_batches_text_output(self, args, capsys):
        service = GCContent(args)
        records = _alignment(
            [
                SeqRecord(Seq("G-C"), id="a"),
                SeqRecord(Seq("A-T"), id="b"),
            ]
        )

        service.calculate_gc_per_sequence(records, is_protein=False)

        out, _ = capsys.readouterr()
        assert out == "a\t1.0\nb\t0.0\n"

    def test_calculate_gc_values_unicode_fallback(self, args):
        service = GCContent(args)
        records = [
            SimpleNamespace(seq="G\u03a9N-", id="a"),
            SimpleNamespace(seq="c?", id="b"),
        ]

        assert service.calculate_gc_per_sequence_data(records, is_protein=False) == [
            ("a", 0.5),
            ("b", 1.0),
        ]
        assert service.calculate_gc_total_value(records, is_protein=False) == round(
            2 / 3,
            4,
        )

    def test_calculate_gc_total_value_mixed_unicode_uses_count_nonzero(
        self, args, mocker
    ):
        service = GCContent(args)
        records = [
            SimpleNamespace(seq="G\u03a9N-", id="a"),
            SimpleNamespace(seq="GC?", id="b"),
        ]
        mocker.patch.object(
            gc_content_module.np,
            "sum",
            side_effect=AssertionError("mixed fallback should use count_nonzero"),
        )

        assert service.calculate_gc_total_value(records, is_protein=False) == round(
            3 / 4,
            4,
        )

    def test_calculate_gc_total_value_mixed_unicode_batches_ascii_fallback(
        self, args, mocker
    ):
        service = GCContent(args)
        records = [
            SimpleNamespace(seq="GCN-", id="a"),
            SimpleNamespace(seq="A-C?", id="b"),
            SimpleNamespace(seq="G\u03a9N-", id="c"),
        ]
        count_nonzero = mocker.spy(gc_content_module.np, "count_nonzero")

        assert service.calculate_gc_total_value(records, is_protein=False) == round(
            4 / 6,
            4,
        )
        assert count_nonzero.call_count == 2

    def test_run_json_summary(self, mocker):
        args = Namespace(fasta="/some/path/to/file.fa", verbose=False, json=True)
        records = _alignment([SeqRecord(Seq("ACGT"), id="a")])
        mocker.patch("phykit.services.alignment.gc_content.get_alignment_and_format", return_value=(records, "fasta", False))
        mocked_json = mocker.patch("phykit.services.alignment.gc_content.print_json")
        service = GCContent(args)
        service.run()
        payload = mocked_json.call_args.args[0]
        assert payload == {"verbose": False, "gc_content": 0.5}

    def test_run_json_verbose(self, mocker):
        args = Namespace(fasta="/some/path/to/file.fa", verbose=True, json=True)
        records = _alignment([SeqRecord(Seq("G-C"), id="a")])
        mocker.patch("phykit.services.alignment.gc_content.get_alignment_and_format", return_value=(records, "fasta", False))
        mocked_json = mocker.patch("phykit.services.alignment.gc_content.print_json")
        service = GCContent(args)
        service.run()
        payload = mocked_json.call_args.args[0]
        assert payload["verbose"] is True
        assert payload["rows"] == [{"taxon": "a", "gc_content": 1.0}]
        assert payload["sequences"] == payload["rows"]

    def test_run_exits_for_protein_alignment(self, mocker):
        args = Namespace(fasta="/some/path/to/file.fa", verbose=False, json=False)
        records = _alignment([SeqRecord(Seq("MSTV"), id="a")])
        mocker.patch("phykit.services.alignment.gc_content.get_alignment_and_format", return_value=(records, "fasta", True))
        mocked_print = mocker.patch("builtins.print")
        service = GCContent(args)
        with pytest.raises(SystemExit) as excinfo:
            service.run()
        assert excinfo.value.code == 2
        mocked_print.assert_called_with("GC content can't be calculated for protein sequences")

    def test_lookup_count_dtype_and_long_invalid_scans(self, monkeypatch):
        gc_content_module._PROTEIN_VALID_LOOKUP = None
        gc_content_module._GC_LOOKUP = None
        protein_lookup = gc_content_module._get_valid_lookup(is_protein=True)
        assert gc_content_module._get_valid_lookup(is_protein=True) is protein_lookup
        gc_lookup = gc_content_module._get_gc_lookup()

        assert not protein_lookup[ord("X")]
        assert protein_lookup[ord("N")]
        assert gc_lookup[ord("G")]
        assert gc_lookup[ord("c")]
        assert not gc_lookup[ord("A")]
        assert gc_content_module._bounded_ascii_count_dtype(
            0x100000000
        ) is gc_content_module.np.uint64

        monkeypatch.setattr(gc_content_module, "_INVALID_SCAN_BYTES", 4)
        assert gc_content_module._has_invalid_bytes(b"-AAAAAAAA", b"-") is True
        assert gc_content_module._has_invalid_bytes(b"AAAAAAAA-", b"-") is True
        assert gc_content_module._has_invalid_bytes(b"AAAA-AAAA", b"-") is True
        assert gc_content_module._has_invalid_bytes(b"AAAAAAAAA", b"-") is False

    def test_ascii_matrix_helper_handles_empty_uneven_and_unicode_records(self):
        assert gc_content_module._gc_counts_from_ascii_matrix(
            [], False
        ) == ([], None, None)

        uneven = [
            SimpleNamespace(seq="ACGT", id="a"),
            SimpleNamespace(seq="ACG", id="b"),
        ]
        record_data, valid_counts, gc_counts = (
            gc_content_module._gc_counts_from_ascii_matrix(uneven, False)
        )
        assert [record_id for record_id, _ in record_data] == ["a", "b"]
        assert valid_counts is None
        assert gc_counts is None

        unicode_records = [
            SimpleNamespace(id="a", seq="AΩGT"),
            SimpleNamespace(id="b", seq="TΩGC"),
        ]
        _, valid_counts, gc_counts = (
            gc_content_module._gc_counts_from_ascii_matrix(
                unicode_records,
                False,
            )
        )
        assert valid_counts is None
        assert gc_counts is None

    def test_forced_per_sequence_fallback_matches_scalar_counts(
        self, args, monkeypatch
    ):
        service = GCContent(args)
        monkeypatch.setattr(
            gc_content_module,
            "_GC_PER_SEQUENCE_SCALAR_MAX_BYTES",
            0,
        )
        records = [
            SimpleNamespace(id="unicode", seq="AΩGC"),
            SimpleNamespace(id="invalid", seq="----"),
        ]

        rows = service.calculate_gc_per_sequence_data(records, is_protein=False)

        assert rows == [("unicode", 0.5), ("invalid", 0.0)]
        monkeypatch.setattr(
            gc_content_module,
            "_GC_PER_SEQUENCE_SCALAR_MAX_BYTES",
            8192,
        )
        assert gc_content_module._gc_per_sequence_data_scalar(
            [SimpleNamespace(id="invalid", seq="----")],
            False,
        ) == [("invalid", 0.0)]

    def test_total_ascii_backends_match_expected_counts(self, monkeypatch):
        assert gc_content_module._gc_total_from_ascii([], False) is None
        records = _alignment(
            [SeqRecord(Seq("ACGT"), id="a"), SeqRecord(Seq("TGCA"), id="b")]
        )
        monkeypatch.setattr(gc_content_module, "_GC_TOTAL_SCALAR_MAX_BYTES", 0)

        assert gc_content_module._gc_total_from_ascii(records, False) == (8, 4)

        invalid_records = _alignment(
            [SeqRecord(Seq("AC-T"), id="a"), SeqRecord(Seq("TG?A"), id="b")]
        )
        monkeypatch.setattr(
            gc_content_module,
            "_GC_TOTAL_BYTE_COUNT_MIN_BYTES",
            1,
        )
        assert gc_content_module._gc_total_from_ascii(
            invalid_records, False
        ) == (6, 2)

        monkeypatch.setattr(
            gc_content_module,
            "_GC_TOTAL_BYTE_COUNT_MIN_BYTES",
            1000,
        )
        assert gc_content_module._gc_total_from_ascii(
            invalid_records, False
        ) == (6, 2)

    def test_run_text_dispatches_verbose_and_summary_paths(self, args, mocker):
        records = _alignment(
            [SeqRecord(Seq("ACGT"), id="a"), SeqRecord(Seq("AGGT"), id="b")]
        )
        service = GCContent(args)
        mocker.patch(
            "phykit.services.alignment.gc_content.get_alignment_and_format",
            return_value=(records, "fasta", False),
        )
        verbose = mocker.patch.object(service, "calculate_gc_per_sequence")
        summary = mocker.patch.object(service, "calculate_gc_total")

        service.verbose = True
        service.run()
        verbose.assert_called_once()
        summary.assert_not_called()

        service.verbose = False
        service.run()
        summary.assert_called_once()

    def test_per_sequence_text_handles_empty_and_broken_pipe(self, args, mocker):
        service = GCContent(args)
        mocker.patch.object(
            service,
            "calculate_gc_per_sequence_data",
            return_value=[],
        )
        mocked_print = mocker.patch("builtins.print")

        service.calculate_gc_per_sequence([], False)
        mocked_print.assert_not_called()

        mocker.patch.object(
            service,
            "calculate_gc_per_sequence_data",
            return_value=[("a", 0.5)],
        )
        mocked_print.side_effect = BrokenPipeError
        service.calculate_gc_per_sequence([], False)

    def test_calculate_gc_total_prints_value(self, args, mocker):
        service = GCContent(args)
        mocked_print = mocker.patch("builtins.print")

        service.calculate_gc_total(
            _alignment([SeqRecord(Seq("ACGT"), id="a")]),
            False,
        )

        mocked_print.assert_called_once_with(0.5)

    def test_total_value_handles_unicode_only_and_empty_records(self, args):
        service = GCContent(args)

        assert service.calculate_gc_total_value(
            [SimpleNamespace(id="unicode", seq="Ω")],
            False,
        ) == 0.0
        with pytest.raises(SystemExit) as exc_info:
            service.calculate_gc_total_value([], False)
        assert exc_info.value.code == 2
