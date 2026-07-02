import pytest
import subprocess
import sys
from argparse import Namespace
from math import isclose
from types import SimpleNamespace

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.occupancy_per_taxon import OccupancyPerTaxon
import phykit.services.alignment.occupancy_per_taxon as occupancy_per_taxon_module


def test_module_import_does_not_import_numpy_or_json_helpers():
    code = """
import sys
import phykit.services.alignment.occupancy_per_taxon as module
assert hasattr(module.np, "__getattr__")
assert callable(module.print_json)
assert "typing" not in sys.modules
assert "numpy" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_lazy_numpy_caches_resolved_attributes():
    lazy_np = occupancy_per_taxon_module._LazyNumpy()

    ones_attr = lazy_np.ones

    assert lazy_np.__dict__["ones"] is ones_attr
    assert lazy_np.ones is ones_attr
    assert lazy_np._module is not None


@pytest.fixture
def args():
    kwargs = dict(alignment="/some/path/to/file.fa")
    return Namespace(**kwargs)


class TestOccupancyPerTaxon(object):
    def test_init_sets_alignment_file_path(self, args):
        occupancy = OccupancyPerTaxon(args)
        assert occupancy.alignment_file_path == args.alignment

    def test_occupancy_per_taxon(self, alignment_simple, args):
        occupancy = OccupancyPerTaxon(args)
        result = occupancy.calculate_occupancy_per_taxon(
            alignment_simple, is_protein=False
        )
        by_taxon = dict(result)

        assert isclose(by_taxon["1"], 5 / 6, rel_tol=0.001)
        assert isclose(by_taxon["2"], 4 / 6, rel_tol=0.001)
        assert isclose(by_taxon["3"], 4 / 6, rel_tol=0.001)
        assert isclose(by_taxon["4"], 5 / 6, rel_tol=0.001)
        assert isclose(by_taxon["5"], 4 / 6, rel_tol=0.001)

    def test_occupancy_per_taxon_handles_lowercase_ambiguity(self, args):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AnXT"), id="a"),
                SeqRecord(Seq("ACxT"), id="b"),
                SeqRecord(Seq("AG-T"), id="c"),
            ]
        )
        occupancy = OccupancyPerTaxon(args)

        dna_result = dict(
            occupancy.calculate_occupancy_per_taxon(alignment, is_protein=False)
        )
        protein_result = dict(
            occupancy.calculate_occupancy_per_taxon(alignment, is_protein=True)
        )

        assert dna_result == {"a": 0.5, "b": 0.75, "c": 0.75}
        assert protein_result == {"a": 0.75, "b": 0.75, "c": 0.75}

    def test_occupancy_per_taxon_uses_ascii_matrix_path(self, args, mocker):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AnXT"), id="a"),
                SeqRecord(Seq("ACxT"), id="b"),
                SeqRecord(Seq("AG-T"), id="c"),
            ]
        )
        occupancy = OccupancyPerTaxon(args)
        mocker.patch.object(
            occupancy_per_taxon_module.np,
            "sum",
            side_effect=AssertionError("ASCII matrix path should avoid row sums"),
        )

        assert dict(
            occupancy.calculate_occupancy_per_taxon(alignment, is_protein=False)
        ) == {"a": 0.5, "b": 0.75, "c": 0.75}

    def test_occupancy_per_taxon_no_gap_ascii_skips_matrix(self, args, mocker):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGT"), id="a"),
                SeqRecord(Seq("TGCA"), id="b"),
            ]
        )
        occupancy = OccupancyPerTaxon(args)
        mocker.patch.object(
            occupancy_per_taxon_module.np,
            "frombuffer",
            side_effect=AssertionError(
                "All-valid ASCII occupancy should avoid matrix construction"
            ),
        )

        assert dict(
            occupancy.calculate_occupancy_per_taxon(alignment, is_protein=False)
        ) == {"a": 1.0, "b": 1.0}

    def test_occupancy_per_taxon_identical_gappy_ascii_skips_matrix(
        self, args, mocker
    ):
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGTN-?X*"), id="a"),
                SeqRecord(Seq("ACGTN-?X*"), id="b"),
                SeqRecord(Seq("ACGTN-?X*"), id="c"),
            ]
        )
        occupancy = OccupancyPerTaxon(args)
        mocker.patch.object(
            occupancy_per_taxon_module.np,
            "frombuffer",
            side_effect=AssertionError(
                "Identical ASCII occupancy should avoid matrix construction"
            ),
        )

        assert dict(
            occupancy.calculate_occupancy_per_taxon(alignment, is_protein=False)
        ) == {"a": 4 / 9, "b": 4 / 9, "c": 4 / 9}
        assert dict(
            occupancy.calculate_occupancy_per_taxon(alignment, is_protein=True)
        ) == {"a": 5 / 9, "b": 5 / 9, "c": 5 / 9}

    def test_occupancy_per_taxon_single_record_skips_matrix_helper(
        self, args, mocker
    ):
        alignment = [SimpleNamespace(seq="ACGTN-?X*", id="solo")]
        occupancy = OccupancyPerTaxon(args)
        mocker.patch(
            "phykit.services.alignment.occupancy_per_taxon._occupancy_from_ascii_matrix",
            side_effect=AssertionError(
                "single-record occupancy should avoid matrix helper setup"
            ),
        )

        assert occupancy.calculate_occupancy_per_taxon(
            alignment,
            is_protein=False,
        ) == [("solo", 4 / 9)]
        assert occupancy.calculate_occupancy_per_taxon(
            alignment,
            is_protein=True,
        ) == [("solo", 5 / 9)]

    def test_occupancy_for_sequence_all_valid_ascii_skips_translate(self):
        class NoTranslateBytes(bytes):
            def translate(self, *_args, **_kwargs):
                raise AssertionError("all-valid ASCII sequence should not translate")

        class EncodedString(str):
            def encode(self, *_args, **_kwargs):
                return NoTranslateBytes(super().encode(*_args, **_kwargs))

        sequence = EncodedString("ACGT" * 2048)

        assert occupancy_per_taxon_module._occupancy_for_sequence(
            sequence,
            is_protein=False,
        ) == 1.0

    def test_invalid_byte_scan_short_buffers_avoid_slice_checks(self):
        class NoSliceBytes(bytes):
            def __getitem__(self, key):
                if isinstance(key, slice):
                    raise AssertionError("short invalid scan should not slice")
                return super().__getitem__(key)

        clean = NoSliceBytes(b"ACGT" * 10)
        gapped = NoSliceBytes(b"ACGTN")

        assert occupancy_per_taxon_module._has_invalid_bytes(
            clean,
            occupancy_per_taxon_module._DNA_INVALID_BYTES,
        ) is False
        assert occupancy_per_taxon_module._has_invalid_bytes(
            gapped,
            occupancy_per_taxon_module._DNA_INVALID_BYTES,
        ) is True

    def test_occupancy_ascii_matrix_identical_rows_return_shared_occupancy(self):
        record_data = [
            ("a", "ACGTN-?X*"),
            ("b", "ACGTN-?X*"),
            ("c", "ACGTN-?X*"),
        ]

        assert occupancy_per_taxon_module._occupancy_from_ascii_matrix(
            record_data,
            is_protein=False,
        ) == [("a", 4 / 9), ("b", 4 / 9), ("c", 4 / 9)]

    def test_occupancy_ascii_matrix_variable_lengths_return_none(self):
        record_data = [
            ("a", "ACGT"),
            ("b", "AC"),
        ]

        assert (
            occupancy_per_taxon_module._occupancy_from_ascii_matrix(
                record_data,
                is_protein=False,
            )
            is None
        )

    def test_occupancy_per_taxon_variable_length_ascii_uses_byte_translate(
        self, args, mocker
    ):
        alignment = [
            SimpleNamespace(seq="ACGTN", id="a"),
            SimpleNamespace(seq="A-C", id="b"),
        ]
        occupancy = OccupancyPerTaxon(args)
        mocker.patch.object(
            occupancy_per_taxon_module.np,
            "frombuffer",
            side_effect=AssertionError(
                "Variable-length ASCII fallback should avoid NumPy per-record setup"
            ),
        )

        assert dict(
            occupancy.calculate_occupancy_per_taxon(alignment, is_protein=False)
        ) == {"a": 0.8, "b": 2 / 3}

    def test_occupancy_per_taxon_unicode_fallback(self, args):
        alignment = [SimpleNamespace(seq="A\u03a9N-", id="a")]
        occupancy = OccupancyPerTaxon(args)

        dna_result = dict(
            occupancy.calculate_occupancy_per_taxon(alignment, is_protein=False)
        )
        protein_result = dict(
            occupancy.calculate_occupancy_per_taxon(alignment, is_protein=True)
        )

        assert dna_result == {"a": 0.5}
        assert protein_result == {"a": 0.75}

    def test_run_text_output_batches_rows(self, mocker, capsys):
        occupancy = OccupancyPerTaxon(Namespace(alignment="x.fa", json=False))
        mocker.patch.object(
            OccupancyPerTaxon,
            "get_alignment_and_format",
            return_value=(object(), "fasta", False),
        )
        mocker.patch.object(
            OccupancyPerTaxon,
            "calculate_occupancy_per_taxon",
            return_value=[("t1", 0.5), ("t2", 1.0)],
        )

        occupancy.run()

        out, _ = capsys.readouterr()
        assert out == "t1\t0.5\nt2\t1.0\n"

    def test_run_json_output_rows(self, mocker):
        occupancy = OccupancyPerTaxon(Namespace(alignment="x.fa", json=True))
        mocker.patch.object(
            OccupancyPerTaxon,
            "get_alignment_and_format",
            return_value=(object(), "fasta", False),
        )
        mocker.patch.object(
            OccupancyPerTaxon,
            "calculate_occupancy_per_taxon",
            return_value=[("t1", 0.55556), ("t2", 1.0)],
        )
        mocked_json = mocker.patch(
            "phykit.services.alignment.occupancy_per_taxon.print_json"
        )

        occupancy.run()

        payload = mocked_json.call_args.args[0]
        assert payload["rows"] == [
            {"taxon": "t1", "occupancy": 0.5556},
            {"taxon": "t2", "occupancy": 1.0},
        ]
        assert payload["taxa"] == payload["rows"]
