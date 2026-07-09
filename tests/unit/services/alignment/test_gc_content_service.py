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

    def test_calculate_gc_per_sequence_data_uses_ascii_matrix_path(self, args, mocker):
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
            side_effect=AssertionError("ASCII matrix path should avoid row sums"),
        )

        assert service.calculate_gc_per_sequence_data(records, is_protein=False) == [
            ("a", 1.0),
            ("b", 0.0),
        ]

    def test_calculate_gc_per_sequence_data_no_gap_ascii_skips_valid_lookup(
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
                "All-valid ASCII GC rows should not build valid lookup counts"
            ),
        )

        assert service.calculate_gc_per_sequence_data(records, is_protein=False) == [
            ("a", 0.75),
            ("b", 0.0),
        ]

    def test_calculate_gc_per_sequence_data_identical_sequences_skip_matrix(
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
                "identical GC rows should not build a byte matrix"
            ),
        )

        assert service.calculate_gc_per_sequence_data(records, is_protein=False) == [
            ("a", 4 / 6),
            ("b", 4 / 6),
            ("c", 4 / 6),
        ]

    def test_calculate_gc_per_sequence_data_variable_length_uses_byte_counts(
        self, args, mocker
    ):
        service = GCContent(args)
        records = [
            SimpleNamespace(seq="GCN-", id="a"),
            SimpleNamespace(seq="A-C", id="b"),
        ]
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
