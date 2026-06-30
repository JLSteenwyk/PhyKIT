from argparse import Namespace
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pytest
import builtins
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.rcvt import RelativeCompositionVariabilityTaxon
import phykit.services.alignment.rcvt as rcvt_module


def _alignment():
    return MultipleSeqAlignment(
        [
            SeqRecord(Seq("ACGT"), id="t1"),
            SeqRecord(Seq("AC-T"), id="t2"),
            SeqRecord(Seq("TCGT"), id="t3"),
        ]
    )


@pytest.fixture
def args():
    return Namespace(alignment="/some/path/to/file.fa")


class TestRCVT:
    def test_module_import_does_not_import_numpy(self):
        code = """
import sys
import phykit.services.alignment.rcvt as module
assert hasattr(module.np, "__getattr__")
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
assert "numpy" not in sys.modules
assert "Bio.Align" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
"""
        subprocess.run([sys.executable, "-c", code], check=True)

    def test_init_sets_expected_attrs(self, args):
        service = RelativeCompositionVariabilityTaxon(args)
        assert service.alignment_file_path == args.alignment
        assert service.json_output is False
        assert service.plot is False
        assert service.plot_output == "rcvt_plot.png"

    def test_calculate_rows_nucleotide(self, args):
        service = RelativeCompositionVariabilityTaxon(args)
        rows = service.calculate_rows(_alignment(), is_protein=False)
        assert len(rows) == 3
        assert rows[0]["taxon"] == "t1"
        assert isinstance(rows[0]["rcvt"], float)

    def test_calculate_rows_preserves_taxon_value_order(self, args):
        service = RelativeCompositionVariabilityTaxon(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AAAAAA"), id="taxon_a"),
                SeqRecord(Seq("AACCCC"), id="taxon_b"),
                SeqRecord(Seq("CCCCCC"), id="taxon_c"),
            ]
        )

        rows = service.calculate_rows(alignment, is_protein=False)

        assert rows == [
            {"taxon": "taxon_a", "rcvt": 0.3704},
            {"taxon": "taxon_b", "rcvt": 0.0741},
            {"taxon": "taxon_c", "rcvt": 0.2963},
        ]

    def test_calculate_rows_empty_alignment(self, args):
        service = RelativeCompositionVariabilityTaxon(args)
        rows = service.calculate_rows(MultipleSeqAlignment([]), is_protein=False)
        assert rows == []

    def test_run_json_with_plot(self, mocker):
        args = Namespace(alignment="/some/path/to/file.fa", json=True, plot=True, plot_output="rcvt.png")
        service = RelativeCompositionVariabilityTaxon(args)
        alignment = _alignment()
        mocker.patch.object(
            RelativeCompositionVariabilityTaxon,
            "get_alignment_and_format",
            return_value=(alignment, "fasta", False),
        )
        mocked_plot = mocker.patch.object(RelativeCompositionVariabilityTaxon, "_plot_rcvt")
        mocked_json = mocker.patch("phykit.services.alignment.rcvt.print_json")

        service.run()

        mocked_plot.assert_called_once()
        payload = mocked_json.call_args.args[0]
        assert payload["rows"] == payload["taxa"]
        assert payload["plot_output"] == "rcvt.png"

    def test_process_args_defaults(self):
        parsed = RelativeCompositionVariabilityTaxon(Namespace(alignment="x.fa")).process_args(
            Namespace(alignment="x.fa")
        )
        assert parsed["json_output"] is False
        assert parsed["plot"] is False
        assert parsed["plot_output"] == "rcvt_plot.png"

    def test_calculate_rows_all_invalid_chars(self, args):
        service = RelativeCompositionVariabilityTaxon(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("----"), id="t1"),
                SeqRecord(Seq("XXXX"), id="t2"),
            ]
        )
        rows = service.calculate_rows(alignment, is_protein=True)
        assert rows == [{"taxon": "t1", "rcvt": 0.0}, {"taxon": "t2", "rcvt": 0.0}]

    def test_calculate_rows_handles_lowercase_ambiguity(self, args):
        service = RelativeCompositionVariabilityTaxon(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AnXT"), id="t1"),
                SeqRecord(Seq("ACxT"), id="t2"),
                SeqRecord(Seq("AG-T"), id="t3"),
            ]
        )

        dna_rows = service.calculate_rows(alignment, is_protein=False)
        protein_rows = service.calculate_rows(alignment, is_protein=True)

        assert [row["rcvt"] for row in dna_rows] == [0.1111, 0.1111, 0.1111]
        assert [row["rcvt"] for row in protein_rows] == [0.1481, 0.1481, 0.1481]

    def test_calculate_rows_ascii_path_uses_lookup(self, mocker, args):
        service = RelativeCompositionVariabilityTaxon(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AnXT"), id="t1"),
                SeqRecord(Seq("ACxT"), id="t2"),
                SeqRecord(Seq("AG-T"), id="t3"),
            ]
        )
        mocker.patch(
            "phykit.services.alignment.rcvt.np.isin",
            side_effect=AssertionError("ASCII path should use the lookup table"),
        )

        rows = service.calculate_rows(alignment, is_protein=False)

        assert [row["rcvt"] for row in rows] == [0.1111, 0.1111, 0.1111]

    def test_calculate_rows_protein_ascii_large_alphabet_uses_bincount(
        self, mocker, args
    ):
        service = RelativeCompositionVariabilityTaxon(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACDEFGHIKLMNPQRSTVWY"), id="t1"),
                SeqRecord(Seq("YWVTSRQPNMLKIHGFEDCA"), id="t2"),
            ]
        )
        bincount_spy = mocker.spy(rcvt_module.np, "bincount")

        rows = service.calculate_rows(alignment, is_protein=True)

        assert bincount_spy.call_count == len(alignment)
        assert [row["taxon"] for row in rows] == ["t1", "t2"]

    def test_calculate_rows_no_gap_ascii_uses_full_lengths(self, mocker, args):
        service = RelativeCompositionVariabilityTaxon(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACDEFGHIKLMNPQRSTVWY"), id="t1"),
                SeqRecord(Seq("ACDEFGHIKLMNPQRSTVWA"), id="t2"),
            ]
        )
        full_spy = mocker.spy(rcvt_module.np, "full")
        mocker.patch(
            "phykit.services.alignment.rcvt.np.isin",
            side_effect=AssertionError("ASCII path should not use Unicode isin"),
        )

        rows = service.calculate_rows(alignment, is_protein=True)

        assert [row["taxon"] for row in rows] == ["t1", "t2"]
        full_spy.assert_called_once_with(
            len(alignment),
            alignment.get_alignment_length(),
            dtype=rcvt_module.np.float64,
        )

    def test_ascii_count_matrix_matches_per_row_bincount_reference(self):
        alphabet = b"ACDEFGHIKLMNPQRSTVWY"
        matrix = rcvt_module.np.frombuffer(
            alphabet + alphabet[::-1] + alphabet[5:] + alphabet[:5],
            dtype=rcvt_module.np.uint8,
        ).reshape(3, len(alphabet))
        unique_chars = rcvt_module.np.unique(matrix)

        observed = RelativeCompositionVariabilityTaxon._ascii_count_matrix(
            matrix,
            unique_chars,
            None,
        )
        expected = rcvt_module.np.zeros(
            (matrix.shape[0], len(unique_chars)),
            dtype=rcvt_module.np.float32,
        )
        for row_idx, row in enumerate(matrix):
            expected[row_idx] = rcvt_module.np.bincount(
                row,
                minlength=256,
            )[unique_chars]

        rcvt_module.np.testing.assert_array_equal(observed, expected)

    def test_ascii_count_matrix_medium_rows_keep_per_row_bincount(self, monkeypatch):
        alphabet = b"ACDEFGHIKLMNPQRSTVWY"
        matrix = rcvt_module.np.tile(
            rcvt_module.np.frombuffer(alphabet, dtype=rcvt_module.np.uint8),
            (8, 20),
        )
        unique_chars = rcvt_module.np.unique(matrix)
        bincount_calls = 0
        original_bincount = rcvt_module.np.bincount

        def count_bincount(*args, **kwargs):
            nonlocal bincount_calls
            bincount_calls += 1
            return original_bincount(*args, **kwargs)

        monkeypatch.setattr(rcvt_module.np, "bincount", count_bincount)

        RelativeCompositionVariabilityTaxon._ascii_count_matrix(
            matrix,
            unique_chars,
            None,
        )

        assert bincount_calls == matrix.shape[0]

    def test_calculate_rows_identical_sequences_skip_matrix_path(self, mocker, args):
        service = RelativeCompositionVariabilityTaxon(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("acgtn-?*"), id="t1"),
                SeqRecord(Seq("ACGTN-?*"), id="t2"),
                SeqRecord(Seq("acgtn-?*"), id="t3"),
            ]
        )
        mocker.patch(
            "phykit.services.alignment.rcvt.np.frombuffer",
            side_effect=AssertionError(
                "identical sequences should not build a NumPy matrix"
            ),
        )

        rows = service.calculate_rows(alignment, is_protein=False)

        assert rows == [
            {"taxon": "t1", "rcvt": 0.0},
            {"taxon": "t2", "rcvt": 0.0},
            {"taxon": "t3", "rcvt": 0.0},
        ]

    def test_identical_sequence_helper_does_not_slice_rows(self):
        class NoSliceList(list):
            def __getitem__(self, key):
                if isinstance(key, slice):
                    raise AssertionError("identical-sequence scan should not slice")
                return super().__getitem__(key)

        sequences = NoSliceList(["ACGT", "ACGT", "ACGT"])

        assert rcvt_module._all_sequences_identical(sequences) is True

    def test_calculate_rows_unicode_fallback(self, args):
        class DummyAlignment(list):
            def get_alignment_length(self):
                return 2

        service = RelativeCompositionVariabilityTaxon(args)
        alignment = DummyAlignment(
            [
                SimpleNamespace(seq="A\u00d1", id="t1"),
                SimpleNamespace(seq="A\u00d1", id="t2"),
                SimpleNamespace(seq="T\u00d1", id="t3"),
            ]
        )

        rows = service.calculate_rows(alignment, is_protein=False)

        assert rows == [
            {"taxon": "t1", "rcvt": 0.1111},
            {"taxon": "t2", "rcvt": 0.1111},
            {"taxon": "t3", "rcvt": 0.2222},
        ]

    def test_plot_rcvt_creates_file(self, tmp_path):
        pytest.importorskip("matplotlib")
        out = tmp_path / "rcvt.png"
        service = RelativeCompositionVariabilityTaxon(
            Namespace(alignment="x.fa", plot=True, plot_output=str(out))
        )
        service._plot_rcvt([{"taxon": "a", "rcvt": 0.1}, {"taxon": "b", "rcvt": 0.2}])
        assert out.exists()

    def test_prepare_rcvt_plot_series_preserves_stable_descending_order(self):
        rows = [
            {"taxon": "a", "rcvt": 0.1},
            {"taxon": "b", "rcvt": 0.3},
            {"taxon": "c", "rcvt": 0.3},
            {"taxon": "d", "rcvt": 0.2},
        ]

        taxa, values = RelativeCompositionVariabilityTaxon._prepare_rcvt_plot_series(
            rows
        )

        assert taxa == ["b", "c", "d", "a"]
        np.testing.assert_allclose(values, np.array([0.3, 0.3, 0.2, 0.1]))

    def test_plot_rcvt_importerror(self, monkeypatch, capsys):
        original_import = builtins.__import__

        def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name.startswith("matplotlib"):
                raise ImportError("no matplotlib")
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", fake_import)
        service = RelativeCompositionVariabilityTaxon(Namespace(alignment="x.fa"))
        with pytest.raises(SystemExit) as exc:
            service._plot_rcvt([{"taxon": "a", "rcvt": 0.1}])
        assert exc.value.code == 2
        out, _ = capsys.readouterr()
        assert "matplotlib is required for --plot in rcvt" in out

    def test_run_text_with_plot(self, mocker, capsys):
        service = RelativeCompositionVariabilityTaxon(
            Namespace(alignment="x.fa", json=False, plot=True, plot_output="rcvt_plot.png")
        )
        mocker.patch.object(
            RelativeCompositionVariabilityTaxon,
            "get_alignment_and_format",
            return_value=(_alignment(), "fasta", False),
        )
        mocker.patch.object(
            RelativeCompositionVariabilityTaxon,
            "calculate_rows",
            return_value=[{"taxon": "a", "rcvt": 0.1}, {"taxon": "b", "rcvt": 0.2}],
        )
        mocked_plot = mocker.patch.object(RelativeCompositionVariabilityTaxon, "_plot_rcvt")
        service.run()
        out, _ = capsys.readouterr()
        assert "a\t0.1\nb\t0.2" in out
        assert "Saved RCVT plot: rcvt_plot.png" in out
        mocked_plot.assert_called_once()
