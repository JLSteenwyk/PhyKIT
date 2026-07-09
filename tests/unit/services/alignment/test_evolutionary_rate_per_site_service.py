from argparse import Namespace
import subprocess
import sys
from types import SimpleNamespace

import pytest
import builtins
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.evolutionary_rate_per_site import EvolutionaryRatePerSite
from phykit.services.alignment.evolutionary_rate_per_site import (
    _ASCII_COUNT_BLOCK_SIZE,
    _GAP_DELETE_TABLES,
    _column_totals,
    _column_sum_squares,
    _column_totals_and_sum_squares_from_ascii_codes,
    _prepare_evolutionary_rate_plot_series,
    _prepare_evolutionary_rate_row_plot_series,
)
import phykit.services.alignment.evolutionary_rate_per_site as erps_module


def _alignment():
    return MultipleSeqAlignment(
        [
            SeqRecord(Seq("ACGT"), id="t1"),
            SeqRecord(Seq("A-GT"), id="t2"),
            SeqRecord(Seq("TCGT"), id="t3"),
        ]
    )


@pytest.fixture
def args():
    return Namespace(alignment="/some/path/to/file.fa")


class TestEvolutionaryRatePerSite:
    def test_module_import_does_not_import_numpy_or_biopython_align(self):
        code = """
import sys
import phykit.services.alignment.evolutionary_rate_per_site as module
assert hasattr(module.np, "__getattr__")
assert callable(module.print_json)
assert "typing" not in sys.modules
assert "numpy" not in sys.modules
assert "Bio.Align" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
"""
        subprocess.run([sys.executable, "-c", code], check=True)

    def test_lazy_numpy_caches_resolved_attributes(self):
        lazy_np = erps_module._LazyNumpy()

        frombuffer_attr = lazy_np.frombuffer

        assert lazy_np.__dict__["frombuffer"] is frombuffer_attr
        assert lazy_np.frombuffer is frombuffer_attr
        assert lazy_np._module is not None

    def test_init_sets_expected_attrs(self, args):
        service = EvolutionaryRatePerSite(args)
        assert service.alignment_file_path == args.alignment
        assert service.json_output is False
        assert service.plot is False
        assert service.plot_output == "evolutionary_rate_per_site_plot.png"

    def test_remove_gap_characters(self, args):
        service = EvolutionaryRatePerSite(args)
        assert service.remove_gap_characters("A-C?T", ["-", "?"]) == "ACT"

    def test_remove_gap_characters_caches_translation_table(self, args):
        service = EvolutionaryRatePerSite(args)
        gap_chars = ["-", "?", "*", "X", "x", "N", "n"]

        _GAP_DELETE_TABLES.pop(tuple(gap_chars), None)
        observed = service.remove_gap_characters("a-c?TnX*", gap_chars)

        assert observed == "ACT"
        assert tuple(gap_chars) in _GAP_DELETE_TABLES

    def test_remove_gap_characters_keeps_multi_character_gap_token_semantics(self, args):
        service = EvolutionaryRatePerSite(args)

        assert service.remove_gap_characters("A--C", ["--"]) == "A--C"

    def test_get_number_of_occurrences_per_character(self, args):
        service = EvolutionaryRatePerSite(args)
        alignment = _alignment()
        counts = service.get_number_of_occurrences_per_character(alignment, 0, ["-"])
        assert counts["A"] == 2
        assert counts["T"] == 1

    def test_get_number_of_occurrences_counts_records_without_column_slice(self, args):
        class NoColumnSliceAlignment:
            def __init__(self, records):
                self.records = records

            def __iter__(self):
                return iter(self.records)

            def __getitem__(self, key):
                if isinstance(key, tuple):
                    raise AssertionError("column slicing should not be used")
                return self.records[key]

        service = EvolutionaryRatePerSite(args)
        alignment = NoColumnSliceAlignment(
            [
                SeqRecord(Seq("A"), id="t1"),
                SeqRecord(Seq("C"), id="t2"),
                SeqRecord(Seq("a"), id="t3"),
                SeqRecord(Seq("-"), id="gap"),
                SeqRecord(Seq("N"), id="ambiguous"),
            ]
        )

        counts = service.get_number_of_occurrences_per_character(
            alignment,
            0,
            ["-", "N", "n"],
        )

        assert counts == {"A": 2, "C": 1}

    def test_get_number_of_occurrences_preserves_multi_character_gap_tokens(self, args):
        service = EvolutionaryRatePerSite(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("A"), id="t1"),
                SeqRecord(Seq("-"), id="dash"),
                SeqRecord(Seq("a"), id="t2"),
            ]
        )

        counts = service.get_number_of_occurrences_per_character(
            alignment,
            0,
            ["--"],
        )

        assert counts == {"A": 2, "-": 1}

    def test_get_number_of_occurrences_multi_character_gaps_avoid_column_slice(
        self, args
    ):
        class NoColumnSliceAlignment:
            def __init__(self, records):
                self.records = records

            def __iter__(self):
                return iter(self.records)

            def __getitem__(self, key):
                if isinstance(key, tuple):
                    raise AssertionError("column slicing should not be used")
                return self.records[key]

        service = EvolutionaryRatePerSite(args)
        alignment = NoColumnSliceAlignment(
            [
                SeqRecord(Seq("A"), id="t1"),
                SeqRecord(Seq("-"), id="dash"),
                SeqRecord(Seq("a"), id="t2"),
            ]
        )

        counts = service.get_number_of_occurrences_per_character(
            alignment,
            0,
            ["--"],
        )

        assert counts == {"A": 2, "-": 1}

    def test_calculate_pic(self, args):
        service = EvolutionaryRatePerSite(args)
        pic = service.calculate_pic({"A": 2, "T": 1})
        assert round(pic, 4) == 0.4444

    def test_calculate_pic_uses_count_sum_of_squares(self, args):
        service = EvolutionaryRatePerSite(args)

        pic = service.calculate_pic({"A": 1, "C": 2, "G": 3, "T": 4})

        assert pic == pytest.approx(0.7)

    def test_calculate_pic_preserves_empty_counts_behavior(self, args):
        service = EvolutionaryRatePerSite(args)

        assert service.calculate_pic({}) == 1

    def test_column_sum_squares_matches_explicit_reduction(self):
        import numpy as np

        counts = np.array(
            [
                [2, 0, 3, 1],
                [1, 4, 0, 2],
                [0, 1, 5, 3],
            ],
            dtype=np.int64,
        )

        np.testing.assert_allclose(
            _column_sum_squares(counts),
            np.sum(counts * counts, axis=0),
        )

    def test_column_totals_use_array_reduction_for_common_counts(self, monkeypatch):
        import numpy as np

        counts = np.arange(4 * 12000, dtype=np.float64).reshape(4, 12000)
        expected = counts.sum(axis=0)

        def fail_sum(*_args, **_kwargs):
            raise AssertionError("common evolutionary-rate totals should use ndarray.sum")

        monkeypatch.setattr(erps_module.np, "sum", fail_sum)

        np.testing.assert_allclose(_column_totals(counts), expected)

    def test_column_totals_preserve_long_matrix_sum_path(self, monkeypatch):
        import numpy as np

        counts = np.ones((4, 20001), dtype=np.float64)
        original_sum = erps_module.np.sum
        calls = []

        def sum_spy(values, *args, **kwargs):
            calls.append((values.shape, kwargs.get("axis")))
            return original_sum(values, *args, **kwargs)

        monkeypatch.setattr(erps_module.np, "sum", sum_spy)

        np.testing.assert_allclose(_column_totals(counts), np.full(20001, 4.0))
        assert calls == [((4, 20001), 0)]

    def test_ascii_column_totals_helper_matches_reference_counts(self):
        import numpy as np

        sequences = [
            "ACDEFGHIKL-?",
            "ACDEFGHIKLXX",
            "CCDEFGHIKLAA",
        ]
        matrix = np.frombuffer(
            "".join(sequences).encode("ascii"),
            dtype=np.uint8,
        ).reshape(len(sequences), len(sequences[0]))
        gap_lookup = erps_module._get_gap_lookup(is_protein=True)
        valid_mask = ~gap_lookup[matrix]
        valid_symbols = np.unique(matrix[valid_mask])

        totals, sum_squares = _column_totals_and_sum_squares_from_ascii_codes(
            matrix,
            valid_mask,
            valid_symbols,
            block_size=5,
        )

        expected_counts = np.array(
            [
                np.count_nonzero(matrix == symbol, axis=0)
                for symbol in valid_symbols
            ],
            dtype=float,
        )
        np.testing.assert_allclose(totals, expected_counts.sum(axis=0))
        np.testing.assert_allclose(sum_squares, (expected_counts ** 2).sum(axis=0))

    def test_ascii_column_totals_helper_counts_all_entries_without_mask(self):
        import numpy as np

        sequences = [
            "ACDEFGHIKL",
            "ACDEFGHIKL",
            "CCDEFGHIKA",
        ]
        matrix = np.frombuffer(
            "".join(sequences).encode("ascii"),
            dtype=np.uint8,
        ).reshape(len(sequences), len(sequences[0]))
        valid_symbols = np.unique(matrix)

        totals, sum_squares = _column_totals_and_sum_squares_from_ascii_codes(
            matrix,
            None,
            valid_symbols,
            block_size=4,
        )

        expected_counts = np.array(
            [
                np.count_nonzero(matrix == symbol, axis=0)
                for symbol in valid_symbols
            ],
            dtype=float,
        )
        np.testing.assert_allclose(totals, expected_counts.sum(axis=0))
        np.testing.assert_allclose(sum_squares, (expected_counts ** 2).sum(axis=0))

    def test_ascii_column_totals_helper_defaults_to_cache_sized_blocks(
        self, monkeypatch
    ):
        import numpy as np

        width = _ASCII_COUNT_BLOCK_SIZE + 1
        matrix = np.tile(
            np.frombuffer(b"ACDEFGHIKLMNPQRSTVWY", dtype=np.uint8),
            (3, width // 20 + 1),
        )[:, :width]
        valid_symbols = np.unique(matrix)
        original_bincount = erps_module.np.bincount
        observed_widths = []

        def bincount_spy(values, *args, **kwargs):
            observed_widths.append(kwargs["minlength"] // 256)
            return original_bincount(values, *args, **kwargs)

        monkeypatch.setattr(erps_module.np, "bincount", bincount_spy)

        _column_totals_and_sum_squares_from_ascii_codes(
            matrix,
            None,
            valid_symbols,
        )

        assert observed_widths == [_ASCII_COUNT_BLOCK_SIZE, 1]

    def test_prepare_evolutionary_rate_plot_series_matches_row_series(self):
        import numpy as np

        values = [0.123456, 0.5, 1.98765]
        rows = [
            {"site": idx + 1, "evolutionary_rate": round(value, 4)}
            for idx, value in enumerate(values)
        ]

        sites, rates = _prepare_evolutionary_rate_plot_series(values)
        row_sites, row_rates = _prepare_evolutionary_rate_row_plot_series(rows)

        np.testing.assert_array_equal(sites, row_sites)
        np.testing.assert_array_equal(rates, row_rates)

    def test_plot_max_uses_array_reduction_for_plot_sized_arrays(self, monkeypatch):
        import numpy as np

        values = np.array([0.1, 0.9, 0.4])

        def fail_max(*_args, **_kwargs):
            raise AssertionError("plot-sized rate maxima should use ndarray.max")

        monkeypatch.setattr(erps_module.np, "max", fail_max)

        assert erps_module._plot_max(values) == pytest.approx(0.9)

    def test_plot_max_preserves_large_array_np_max_path(self, monkeypatch):
        import numpy as np

        values = np.ones(erps_module._PLOT_DIRECT_MAX_LIMIT + 1, dtype=np.float64)
        original_max = erps_module.np.max
        calls = []

        def max_spy(observed, *args, **kwargs):
            calls.append(observed.shape)
            return original_max(observed, *args, **kwargs)

        monkeypatch.setattr(erps_module.np, "max", max_spy)

        assert erps_module._plot_max(values) == pytest.approx(1.0)
        assert calls == [(erps_module._PLOT_DIRECT_MAX_LIMIT + 1,)]

    def test_calculate_evolutionary_rate_per_site(self, args):
        service = EvolutionaryRatePerSite(args)
        values = service.calculate_evolutionary_rate_per_site(_alignment(), is_protein=False)
        assert len(values) == 4
        assert all(0 <= v <= 1 for v in values)

    def test_run_json_with_plot(self, mocker):
        args = Namespace(alignment="/some/path/to/file.fa", json=True, plot=True, plot_output="erps.png")
        service = EvolutionaryRatePerSite(args)
        alignment = _alignment()
        mocker.patch.object(
            EvolutionaryRatePerSite,
            "get_alignment_and_format",
            return_value=(alignment, "fasta", False),
        )
        mocked_plot = mocker.patch.object(EvolutionaryRatePerSite, "_plot_evolutionary_rate_per_site")
        mocked_json = mocker.patch("phykit.services.alignment.evolutionary_rate_per_site.print_json")

        service.run()

        mocked_plot.assert_called_once()
        payload = mocked_json.call_args.args[0]
        assert payload["rows"] == payload["sites"]
        assert payload["plot_output"] == "erps.png"

    def test_run_json_without_plot_builds_site_rows(self, mocker):
        service = EvolutionaryRatePerSite(
            Namespace(alignment="x.fa", json=True, plot=False)
        )
        mocker.patch.object(
            EvolutionaryRatePerSite,
            "get_alignment_and_format",
            return_value=("alignment", "fasta", False),
        )
        mocker.patch.object(
            EvolutionaryRatePerSite,
            "calculate_evolutionary_rate_per_site",
            return_value=[0.0, 1.23456],
        )
        mocked_json = mocker.patch(
            "phykit.services.alignment.evolutionary_rate_per_site.print_json"
        )

        service.run()

        payload = mocked_json.call_args.args[0]
        assert payload == {
            "rows": [
                {"site": 1, "evolutionary_rate": 0.0},
                {"site": 2, "evolutionary_rate": 1.2346},
            ],
            "sites": [
                {"site": 1, "evolutionary_rate": 0.0},
                {"site": 2, "evolutionary_rate": 1.2346},
            ],
        }

    def test_run_text_without_plot_formats_site_values_directly(self, mocker, capsys):
        service = EvolutionaryRatePerSite(Namespace(alignment="x.fa", json=False, plot=False))
        mocker.patch.object(
            EvolutionaryRatePerSite,
            "get_alignment_and_format",
            return_value=("alignment", "fasta", False),
        )
        mocker.patch.object(
            EvolutionaryRatePerSite,
            "calculate_evolutionary_rate_per_site",
            return_value=[0.0, 1.23456],
        )

        service.run()

        out, _ = capsys.readouterr()
        assert out == "1\t0.0\n2\t1.2346\n"

    def test_process_args_defaults(self):
        parsed = EvolutionaryRatePerSite(Namespace(alignment="x.fa")).process_args(
            Namespace(alignment="x.fa")
        )
        assert parsed["json_output"] is False
        assert parsed["plot"] is False
        assert parsed["plot_output"] == "evolutionary_rate_per_site_plot.png"

    def test_calculate_evolutionary_rate_per_site_all_gaps(self, args):
        service = EvolutionaryRatePerSite(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("----"), id="t1"),
                SeqRecord(Seq("????"), id="t2"),
            ]
        )
        values = service.calculate_evolutionary_rate_per_site(alignment, is_protein=False)
        assert values == [0, 0, 0, 0]

    def test_calculate_evolutionary_rate_per_site_single_valid_symbol_shortcut(
        self, mocker, args
    ):
        service = EvolutionaryRatePerSite(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AAAA"), id="t1"),
                SeqRecord(Seq("A-A?"), id="t2"),
                SeqRecord(Seq("AAAA"), id="t3"),
            ]
        )
        mocker.patch(
            "phykit.services.alignment.evolutionary_rate_per_site.np.array",
            side_effect=AssertionError("single-symbol PIC should skip counts"),
        )

        values = service.calculate_evolutionary_rate_per_site(
            alignment,
            is_protein=False,
        )

        assert values == [0.0, 0.0, 0.0, 0.0]

    def test_calculate_evolutionary_rate_per_site_identical_sequences_skip_matrix(
        self, mocker, args
    ):
        service = EvolutionaryRatePerSite(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AcGt"), id="t1"),
                SeqRecord(Seq("aCgT"), id="t2"),
                SeqRecord(Seq("ACGT"), id="t3"),
            ]
        )
        mocker.patch(
            "phykit.services.alignment.evolutionary_rate_per_site.np.frombuffer",
            side_effect=AssertionError("identical alignments should skip matrix setup"),
        )

        values = service.calculate_evolutionary_rate_per_site(
            alignment,
            is_protein=False,
        )

        assert values == [0.0, 0.0, 0.0, 0.0]

    def test_identical_sequence_helper_does_not_slice_rows(self):
        class NoSliceList(list):
            def __getitem__(self, key):
                if isinstance(key, slice):
                    raise AssertionError("identical-sequence scan should not slice")
                return super().__getitem__(key)

        sequences = NoSliceList(["ACGT", "ACGT", "ACGT"])

        assert erps_module._all_sequences_identical(sequences) is True

    def test_calculate_evolutionary_rate_per_site_handles_lowercase_ambiguity(self, args):
        service = EvolutionaryRatePerSite(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AnXT"), id="t1"),
                SeqRecord(Seq("ACxT"), id="t2"),
                SeqRecord(Seq("AG-T"), id="t3"),
            ]
        )

        dna_values = service.calculate_evolutionary_rate_per_site(alignment, is_protein=False)
        protein_values = service.calculate_evolutionary_rate_per_site(alignment, is_protein=True)

        assert dna_values == [0.0, 0.5, 0.0, 0.0]
        assert protein_values[0] == 0.0
        assert round(protein_values[1], 4) == 0.6667
        assert protein_values[2:] == [0.0, 0.0]

    def test_calculate_evolutionary_rate_per_site_ascii_path_uses_lookup(self, mocker, args):
        service = EvolutionaryRatePerSite(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AnXT"), id="t1"),
                SeqRecord(Seq("ACxT"), id="t2"),
                SeqRecord(Seq("AG-T"), id="t3"),
            ]
        )
        mocker.patch(
            "phykit.services.alignment.evolutionary_rate_per_site.np.isin",
            side_effect=AssertionError("ASCII path should use the lookup table"),
        )

        values = service.calculate_evolutionary_rate_per_site(
            alignment,
            is_protein=False,
        )

        assert values == [0.0, 0.5, 0.0, 0.0]

    def test_calculate_evolutionary_rate_per_site_small_alphabet_uses_array_sum(
        self,
        monkeypatch,
        args,
    ):
        service = EvolutionaryRatePerSite(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGT"), id="t1"),
                SeqRecord(Seq("ACGA"), id="t2"),
                SeqRecord(Seq("TCGT"), id="t3"),
            ]
        )

        def fail_sum(*_args, **_kwargs):
            raise AssertionError("small-alphabet ERPS counts should use ndarray.sum")

        monkeypatch.setattr(erps_module.np, "sum", fail_sum)

        values = service.calculate_evolutionary_rate_per_site(
            alignment,
            is_protein=False,
        )

        assert values == pytest.approx([4 / 9, 0.0, 0.0, 4 / 9])

    def test_calculate_evolutionary_rate_per_site_protein_uses_block_counts(
        self, mocker, args
    ):
        service = EvolutionaryRatePerSite(args)
        alphabet = "ACDEFGHIKL"
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq(symbol + "A"), id=f"t{idx}")
                for idx, symbol in enumerate(alphabet)
            ]
        )
        mocker.patch(
            "phykit.services.alignment.evolutionary_rate_per_site.np.isin",
            side_effect=AssertionError("ASCII protein path should use byte counts"),
        )
        spy = mocker.spy(erps_module, "_column_totals_and_sum_squares_from_ascii_codes")

        values = service.calculate_evolutionary_rate_per_site(
            alignment,
            is_protein=True,
        )

        assert spy.call_args.args[1] is None
        assert values[0] == pytest.approx(0.9)
        assert values[1] == 0.0

    def test_calculate_evolutionary_rate_per_site_unicode_fallback(self, args):
        class DummyAlignment(list):
            def get_alignment_length(self):
                return 2

        service = EvolutionaryRatePerSite(args)
        alignment = DummyAlignment(
            [
                SimpleNamespace(seq="A\u00d1", id="t1"),
                SimpleNamespace(seq="A\u00d1", id="t2"),
                SimpleNamespace(seq="T\u00d1", id="t3"),
            ]
        )

        values = service.calculate_evolutionary_rate_per_site(
            alignment,
            is_protein=False,
        )

        assert round(values[0], 4) == 0.4444
        assert values[1] == 0.0

    def test_plot_evolutionary_rate_per_site_creates_file(self, tmp_path):
        pytest.importorskip("matplotlib")
        out = tmp_path / "erps.png"
        service = EvolutionaryRatePerSite(
            Namespace(alignment="x.fa", plot=True, plot_output=str(out))
        )
        service._plot_evolutionary_rate_per_site(
            [{"site": 1, "evolutionary_rate": 0.2}, {"site": 2, "evolutionary_rate": 0.8}]
        )
        assert out.exists()

    def test_plot_evolutionary_rate_skips_redundant_tight_layout(
        self, tmp_path, monkeypatch
    ):
        pytest.importorskip("matplotlib")
        from matplotlib.figure import Figure

        out = tmp_path / "erps.png"
        service = EvolutionaryRatePerSite(
            Namespace(alignment="x.fa", plot=True, plot_output=str(out))
        )

        def fail_tight_layout(self, *args, **kwargs):
            raise AssertionError("bbox_inches='tight' handles saved bounds")

        monkeypatch.setattr(Figure, "tight_layout", fail_tight_layout)
        service._plot_evolutionary_rate_per_site(
            [
                {"site": 1, "evolutionary_rate": 0.2},
                {"site": 2, "evolutionary_rate": 0.4},
            ]
        )

        assert out.exists()

    def test_plot_evolutionary_rate_per_site_empty_rows(self):
        pytest.importorskip("matplotlib")
        service = EvolutionaryRatePerSite(Namespace(alignment="x.fa", plot=True))
        service._plot_evolutionary_rate_per_site([])

    def test_plot_evolutionary_rate_per_site_importerror(self, monkeypatch, capsys):
        original_import = builtins.__import__

        def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name.startswith("matplotlib"):
                raise ImportError("no matplotlib")
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", fake_import)
        service = EvolutionaryRatePerSite(Namespace(alignment="x.fa"))
        with pytest.raises(SystemExit) as exc:
            service._plot_evolutionary_rate_per_site([{"site": 1, "evolutionary_rate": 0.2}])
        assert exc.value.code == 2
        out, _ = capsys.readouterr()
        assert "matplotlib is required for --plot in evolutionary_rate_per_site" in out

    def test_run_text_with_plot_message(self, mocker, capsys):
        service = EvolutionaryRatePerSite(
            Namespace(alignment="x.fa", json=False, plot=True, plot_output="erps.png")
        )
        mocker.patch.object(
            EvolutionaryRatePerSite,
            "get_alignment_and_format",
            return_value=(_alignment(), "fasta", False),
        )
        mocker.patch.object(
            EvolutionaryRatePerSite,
            "calculate_evolutionary_rate_per_site",
            return_value=[0.2, 0.4],
        )
        mocked_plot = mocker.patch.object(EvolutionaryRatePerSite, "_plot_evolutionary_rate_values")
        service.run()
        out, _ = capsys.readouterr()
        assert "1\t0.2\n2\t0.4" in out
        assert "Saved evolutionary-rate plot: erps.png" in out
        mocked_plot.assert_called_once_with([0.2, 0.4])
