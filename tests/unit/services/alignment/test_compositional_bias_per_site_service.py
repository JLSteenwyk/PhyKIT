from argparse import Namespace
import subprocess
import sys

import numpy as np
import pytest
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.compositional_bias_per_site import (
    CompositionalBiasPerSite,
)
import phykit.services.alignment.compositional_bias_per_site as cbps_module
import builtins


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


class TestCompositionalBiasPerSite:
    def test_module_import_does_not_import_numpy_or_biopython_align(self):
        code = """
import sys
import phykit.services.alignment.compositional_bias_per_site as module
assert hasattr(module.np, "__getattr__")
assert "typing" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
assert "numpy" not in sys.modules
assert "Bio.Align" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
"""
        subprocess.run([sys.executable, "-c", code], check=True)

    def test_lazy_numpy_caches_resolved_attributes(self):
        lazy_np = cbps_module._LazyNumpy()

        array_attr = lazy_np.array

        assert lazy_np.__dict__["array"] is array_attr
        assert lazy_np.array is array_attr
        assert lazy_np._module is not None

    def test_init_sets_expected_attrs(self, args):
        service = CompositionalBiasPerSite(args)
        assert service.alignment_file_path == args.alignment
        assert service.json_output is False
        assert service.plot is False
        assert service.plot_output == "compositional_bias_per_site_plot.png"

    def test_get_number_of_occurrences_per_character(self, args):
        service = CompositionalBiasPerSite(args)
        counts = service.get_number_of_occurrences_per_character(_alignment(), 0, is_protein=False)
        assert sorted(counts) == [1, 2]

    def test_get_number_of_occurrences_preserves_first_seen_count_order(self, args):
        service = CompositionalBiasPerSite(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("A"), id="t1"),
                SeqRecord(Seq("C"), id="t2"),
                SeqRecord(Seq("A"), id="t3"),
                SeqRecord(Seq("G"), id="t4"),
                SeqRecord(Seq("C"), id="t5"),
                SeqRecord(Seq("C"), id="t6"),
                SeqRecord(Seq("-"), id="gap"),
            ]
        )

        counts = service.get_number_of_occurrences_per_character(
            alignment, 0, is_protein=False
        )

        assert counts == [2, 3, 1]

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

        service = CompositionalBiasPerSite(args)
        alignment = NoColumnSliceAlignment(
            [
                SeqRecord(Seq("A"), id="t1"),
                SeqRecord(Seq("C"), id="t2"),
                SeqRecord(Seq("A"), id="t3"),
                SeqRecord(Seq("G"), id="t4"),
                SeqRecord(Seq("C"), id="t5"),
                SeqRecord(Seq("C"), id="t6"),
                SeqRecord(Seq("-"), id="gap"),
            ]
        )

        counts = service.get_number_of_occurrences_per_character(
            alignment, 0, is_protein=False
        )

        assert counts == [2, 3, 1]

    def test_calculate_compositional_bias_per_site(self, args):
        service = CompositionalBiasPerSite(args)
        stat_res, corrected = service.calculate_compositional_bias_per_site(_alignment(), is_protein=False)
        assert len(stat_res) == 4
        assert len(corrected) == 4

    def test_false_discovery_control_matches_bh_adjustment(self):
        corrected = cbps_module._false_discovery_control([0.01, 0.04, 0.03, 0.2])
        assert corrected == pytest.approx([0.04, 0.0533333333, 0.0533333333, 0.2])

    def test_medium_false_discovery_control_matches_scalar_reference(self):
        p_values = [((idx * 37) % 101) / 1000 for idx in range(32)]
        indexed = sorted(enumerate(p_values), key=lambda x: x[1])
        expected = [0.0] * len(p_values)
        previous = 1.0
        for rank_minus_1 in range(len(p_values) - 1, -1, -1):
            original_idx, p_value = indexed[rank_minus_1]
            rank = rank_minus_1 + 1
            adjusted = min(p_value * len(p_values) / rank, previous)
            adjusted = min(adjusted, 1.0)
            expected[original_idx] = adjusted
            previous = adjusted

        corrected = cbps_module._false_discovery_control(p_values)
        assert corrected == pytest.approx(expected)

    def test_vector_false_discovery_control_caps_adjusted_values(self):
        p_values = [1.0 if idx % 5 else 0.001 * (idx + 1) for idx in range(64)]
        indexed = sorted(enumerate(p_values), key=lambda x: x[1])
        expected = [0.0] * len(p_values)
        previous = 1.0
        for rank_minus_1 in range(len(p_values) - 1, -1, -1):
            original_idx, p_value = indexed[rank_minus_1]
            rank = rank_minus_1 + 1
            adjusted = min(p_value * len(p_values) / rank, previous)
            adjusted = min(adjusted, 1.0)
            expected[original_idx] = adjusted
            previous = adjusted

        corrected = cbps_module._false_discovery_control(p_values)
        assert corrected == pytest.approx(expected)
        assert max(expected) == 1.0

    def test_erfc_array_vectorizes_square_root(self, mocker):
        values = np.array([0.0, 0.5, 2.0])
        sqrt_spy = mocker.spy(cbps_module.np, "sqrt")

        observed = cbps_module._erfc_array(values)

        sqrt_spy.assert_called_once()
        np.testing.assert_allclose(sqrt_spy.call_args.args[0], values)
        expected = np.array([cbps_module.erfc(np.sqrt(value)) for value in values])
        np.testing.assert_allclose(observed, expected)

    def test_small_false_discovery_control_does_not_import_numpy(self):
        code = """
import sys
import phykit.services.alignment.compositional_bias_per_site as module
corrected = module._false_discovery_control([0.01, 0.04, 0.03, 0.2])
assert [round(value, 10) for value in corrected] == [0.04, 0.0533333333, 0.0533333333, 0.2]
assert "numpy" not in sys.modules
"""
        subprocess.run([sys.executable, "-c", code], check=True)

    def test_calculate_compositional_bias_per_site_does_not_import_scipy(self, monkeypatch, args):
        original_import = builtins.__import__

        def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name == "scipy" or name.startswith("scipy."):
                raise AssertionError("compositional bias calculation should not import scipy")
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", fake_import)
        service = CompositionalBiasPerSite(args)

        stat_res, corrected = service.calculate_compositional_bias_per_site(
            _alignment(),
            is_protein=False,
        )

        assert round(float(stat_res[0].pvalue), 4) == 0.5637
        assert corrected[0] == pytest.approx(float(stat_res[0].pvalue))

    def test_build_rows_handles_nan(self, args):
        service = CompositionalBiasPerSite(args)
        stat_res = [
            type("R", (), {"statistic": 0.0, "pvalue": np.nan})(),
            type("R", (), {"statistic": 0.2, "pvalue": 0.6547})(),
        ]
        rows = service._build_rows(stat_res, ["nan", 1.0])
        assert rows[0] == {
            "site": 1,
            "chi_square": 0.0,
            "p_value_corrected": None,
            "p_value": None,
        }
        assert rows[1] == {
            "site": 2,
            "chi_square": 0.2,
            "p_value_corrected": 1.0,
            "p_value": 0.6547,
        }

    def test_prepare_compositional_bias_plot_series_matches_row_series(self, args):
        service = CompositionalBiasPerSite(args)
        stat_res = [
            type("R", (), {"statistic": 0.0, "pvalue": np.nan})(),
            type("R", (), {"statistic": 0.2, "pvalue": 0.6547})(),
            type("R", (), {"statistic": 0.4, "pvalue": 0.00123})(),
        ]
        corrected = ["nan", 0.987654, 0.00123]
        rows = service._build_rows(stat_res, corrected)

        sites, pvals = cbps_module._prepare_compositional_bias_plot_series(corrected)
        row_sites, row_pvals = cbps_module._prepare_compositional_bias_row_plot_series(rows)

        np.testing.assert_array_equal(sites, row_sites)
        np.testing.assert_allclose(pvals, row_pvals, equal_nan=True)

    def test_plot_max_uses_array_reduction_for_plot_sized_arrays(self, monkeypatch):
        values = np.array([1, 3, 2], dtype=np.int32)

        def fail_max(*_args, **_kwargs):
            raise AssertionError("plot-sized compositional-bias maxima should use ndarray.max")

        monkeypatch.setattr(cbps_module.np, "max", fail_max)

        assert cbps_module._plot_max(values) == 3

    def test_plot_max_preserves_large_array_np_max_path(self, monkeypatch):
        values = np.ones(cbps_module._PLOT_DIRECT_MAX_LIMIT + 1, dtype=np.int32)
        original_max = cbps_module.np.max
        calls = []

        def max_spy(observed, *args, **kwargs):
            calls.append(observed.shape)
            return original_max(observed, *args, **kwargs)

        monkeypatch.setattr(cbps_module.np, "max", max_spy)

        assert cbps_module._plot_max(values) == 1
        assert calls == [(cbps_module._PLOT_DIRECT_MAX_LIMIT + 1,)]

    def test_power_divergence_results_from_arrays_preserves_float_results(self):
        statistics = np.array([0, 1.25, 2.5], dtype=np.float64)
        p_values = np.array([np.nan, 0.5, 0.125], dtype=np.float64)

        results = cbps_module._power_divergence_results_from_arrays(
            statistics,
            p_values,
        )

        assert [type(result) for result in results] == [
            cbps_module.Power_divergenceResult,
            cbps_module.Power_divergenceResult,
            cbps_module.Power_divergenceResult,
        ]
        assert [result.statistic for result in results] == [0.0, 1.25, 2.5]
        assert np.isnan(results[0].pvalue)
        assert [result.pvalue for result in results[1:]] == [0.5, 0.125]

    def test_pvalues_for_fdr_all_valid_uses_direct_tolist(self):
        class PValues:
            def __getitem__(self, key):
                raise AssertionError("all-valid p-values should not be masked")

            def tolist(self):
                return [0.1, 0.2, 0.3]

        class NoNanMask:
            def any(self):
                return False

        assert cbps_module._pvalues_for_fdr(PValues(), NoNanMask()) == [
            0.1,
            0.2,
            0.3,
        ]

    def test_pvalues_for_fdr_mixed_nan_filters_values(self):
        p_values = np.array([np.nan, 0.5, 0.125], dtype=np.float64)
        nan_mask = np.isnan(p_values)

        assert cbps_module._pvalues_for_fdr(p_values, nan_mask) == [0.5, 0.125]

    def test_restore_nan_corrected_pvalues_all_valid_returns_directly(self, mocker):
        corrected = [0.2, 0.4, 0.6]
        nan_mask = np.array([False, False, False], dtype=bool)
        mocker.patch.object(
            cbps_module.np,
            "flatnonzero",
            side_effect=AssertionError("all-valid p-values should not rebuild slots"),
        )

        restored = cbps_module._restore_nan_corrected_pvalues(corrected, nan_mask)

        assert restored is corrected

    def test_run_json_with_plot(self, mocker):
        args = Namespace(alignment="/some/path/to/file.fa", json=True, plot=True, plot_output="cbps.png")
        service = CompositionalBiasPerSite(args)
        mocker.patch.object(
            CompositionalBiasPerSite,
            "get_alignment_and_format",
            return_value=(_alignment(), "fasta", False),
        )
        mocker.patch.object(
            CompositionalBiasPerSite,
            "calculate_compositional_bias_per_site",
            return_value=(
                [type("R", (), {"statistic": 0.2, "pvalue": 0.6547})()],
                [1.0],
            ),
        )
        mocked_plot = mocker.patch.object(CompositionalBiasPerSite, "_plot_compositional_bias_manhattan")
        mocked_json = mocker.patch("phykit.services.alignment.compositional_bias_per_site.print_json")

        service.run()

        mocked_plot.assert_called_once()
        payload = mocked_json.call_args.args[0]
        assert payload["rows"] == payload["sites"]
        assert payload["plot_output"] == "cbps.png"

    def test_run_text_without_plot_formats_rows_directly(self, mocker, capsys):
        service = CompositionalBiasPerSite(
            Namespace(alignment="x.fa", json=False, plot=False)
        )
        mocker.patch.object(
            CompositionalBiasPerSite,
            "get_alignment_and_format",
            return_value=("alignment", "fasta", False),
        )
        mocker.patch.object(
            CompositionalBiasPerSite,
            "calculate_compositional_bias_per_site",
            return_value=(
                [
                    type("R", (), {"statistic": 2.0, "pvalue": 0.0456})(),
                    type("R", (), {"statistic": 0.0, "pvalue": np.nan})(),
                ],
                [0.1234, "nan"],
            ),
        )
        mocker.patch.object(
            service,
            "_build_rows",
            side_effect=AssertionError("text output should not build row dictionaries"),
        )

        service.run()

        out, _ = capsys.readouterr()
        assert out == "1\t2.0\t0.1234\t0.0456\n2\t0.0\tnan\tnan\n"

    def test_format_text_rows_preserves_order_rounding_and_nan_values(self):
        rows = CompositionalBiasPerSite._format_text_rows(
            [
                cbps_module.Power_divergenceResult(2.12345, 0.04567),
                cbps_module.Power_divergenceResult(0.0, np.nan),
            ],
            [0.12345, "nan"],
        )

        assert rows == "1\t2.1235\t0.1235\t0.0457\n2\t0.0\tnan\tnan"

    def test_process_args_defaults(self):
        parsed = CompositionalBiasPerSite(Namespace(alignment="x.fa")).process_args(
            Namespace(alignment="x.fa")
        )
        assert parsed["json_output"] is False
        assert parsed["plot"] is False
        assert parsed["plot_output"] == "compositional_bias_per_site_plot.png"
        assert parsed["plot_config"] is None

    def test_no_plot_init_does_not_import_plot_config(self):
        code = """
import sys
from argparse import Namespace
from phykit.services.alignment.compositional_bias_per_site import CompositionalBiasPerSite
CompositionalBiasPerSite(Namespace(alignment="x.fa", plot=False))
assert "phykit.helpers.plot_config" not in sys.modules
"""
        subprocess.run([sys.executable, "-c", code], check=True)

    def test_process_args_builds_plot_config_when_plotting(self):
        parsed = CompositionalBiasPerSite(
            Namespace(alignment="x.fa", plot=True)
        ).process_args(
            Namespace(alignment="x.fa", plot=True)
        )

        assert parsed["plot"] is True
        assert parsed["plot_config"] is not None

    def test_calculate_compositional_bias_per_site_all_gaps(self):
        service = CompositionalBiasPerSite(Namespace(alignment="x.fa"))
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("----"), id="t1"),
                SeqRecord(Seq("????"), id="t2"),
            ]
        )
        stat_res, corrected = service.calculate_compositional_bias_per_site(alignment, is_protein=False)
        assert len(stat_res) == 4
        assert corrected == ["nan", "nan", "nan", "nan"]

    def test_calculate_compositional_bias_per_site_single_valid_symbol_shortcut(
        self, mocker
    ):
        service = CompositionalBiasPerSite(Namespace(alignment="x.fa"))
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AAAA"), id="t1"),
                SeqRecord(Seq("A-A?"), id="t2"),
                SeqRecord(Seq("AAAA"), id="t3"),
            ]
        )
        mocker.patch(
            "phykit.services.alignment.compositional_bias_per_site.np.array",
            side_effect=AssertionError("single-symbol bias should skip counts"),
        )

        stat_res, corrected = service.calculate_compositional_bias_per_site(
            alignment,
            is_protein=False,
        )

        assert [result.statistic for result in stat_res] == [0.0, 0.0, 0.0, 0.0]
        assert all(np.isnan(result.pvalue) for result in stat_res)
        assert corrected == ["nan", "nan", "nan", "nan"]

    def test_calculate_compositional_bias_per_site_identical_sequences_skip_matrix(
        self, mocker
    ):
        service = CompositionalBiasPerSite(Namespace(alignment="x.fa"))
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AcGt"), id="t1"),
                SeqRecord(Seq("aCgT"), id="t2"),
                SeqRecord(Seq("ACGT"), id="t3"),
                SeqRecord(Seq("acgt"), id="t4"),
            ]
        )
        mocker.patch(
            "phykit.services.alignment.compositional_bias_per_site.np.frombuffer",
            side_effect=AssertionError("identical alignments should skip matrix setup"),
        )

        stat_res, corrected = service.calculate_compositional_bias_per_site(
            alignment,
            is_protein=False,
        )

        assert [result.statistic for result in stat_res] == [0.0, 0.0, 0.0, 0.0]
        assert all(np.isnan(result.pvalue) for result in stat_res)
        assert corrected == ["nan", "nan", "nan", "nan"]

    def test_calculate_compositional_bias_per_site_small_alignment_uses_scalar_path(
        self, mocker, args
    ):
        service = CompositionalBiasPerSite(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGTAC"), id="t1"),
                SeqRecord(Seq("A-GTCC"), id="t2"),
                SeqRecord(Seq("TCGTAG"), id="t3"),
            ]
        )
        mocker.patch(
            "phykit.services.alignment.compositional_bias_per_site.np.frombuffer",
            side_effect=AssertionError("small alignments should use scalar counts"),
        )

        stat_res, corrected = service.calculate_compositional_bias_per_site(
            alignment,
            is_protein=False,
        )

        assert [round(float(result.statistic), 4) for result in stat_res] == [
            0.3333,
            0.0,
            0.0,
            0.0,
            0.3333,
            0.3333,
        ]
        assert corrected[1] == "nan"
        assert all(
            value == pytest.approx(0.5637028616507729)
            for value in corrected
            if value != "nan"
        )

    def test_scalar_compositional_bias_matches_vector_path_for_small_input(self, args):
        service = CompositionalBiasPerSite(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGTAC"), id="t1"),
                SeqRecord(Seq("A-GTCC"), id="t2"),
                SeqRecord(Seq("TCGTAG"), id="t3"),
            ]
        )
        raw_sequences = [str(record.seq) for record in alignment]
        scalar_stat, scalar_corrected = cbps_module._calculate_compositional_bias_scalar(
            raw_sequences,
            alignment.get_alignment_length(),
            is_protein=False,
        )
        original_limit = cbps_module._SCALAR_CELL_MAX
        cbps_module._SCALAR_CELL_MAX = 0
        try:
            vector_stat, vector_corrected = service.calculate_compositional_bias_per_site(
                alignment,
                is_protein=False,
            )
        finally:
            cbps_module._SCALAR_CELL_MAX = original_limit

        assert [result.statistic for result in scalar_stat] == pytest.approx(
            [result.statistic for result in vector_stat]
        )
        assert [result.pvalue for result in scalar_stat] == pytest.approx(
            [result.pvalue for result in vector_stat],
            nan_ok=True,
        )
        for scalar_value, vector_value in zip(scalar_corrected, vector_corrected):
            if scalar_value == "nan" or vector_value == "nan":
                assert scalar_value == vector_value
            else:
                assert scalar_value == pytest.approx(vector_value)

    def test_corrected_p_values_keep_interleaved_nan_positions(self, args):
        service = CompositionalBiasPerSite(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AA-ACN"), id="t1"),
                SeqRecord(Seq("AC-CCN"), id="t2"),
                SeqRecord(Seq("AA?GAN"), id="t3"),
            ]
        )

        _stat_res, corrected = service.calculate_compositional_bias_per_site(
            alignment,
            is_protein=False,
        )

        assert corrected[0] == "nan"
        assert corrected[2] == "nan"
        assert corrected[5] == "nan"
        assert all(isinstance(corrected[idx], float) for idx in (1, 3, 4))

    def test_calculate_compositional_bias_per_site_handles_lowercase_ambiguity(self, args):
        service = CompositionalBiasPerSite(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AnXT"), id="t1"),
                SeqRecord(Seq("ACxT"), id="t2"),
                SeqRecord(Seq("AG-T"), id="t3"),
            ]
        )

        dna_stat_res, dna_corrected = service.calculate_compositional_bias_per_site(
            alignment,
            is_protein=False,
        )
        protein_stat_res, protein_corrected = service.calculate_compositional_bias_per_site(
            alignment,
            is_protein=True,
        )

        assert [float(result.statistic) for result in dna_stat_res] == [0.0, 0.0, 0.0, 0.0]
        assert dna_corrected == ["nan", 1.0, "nan", "nan"]
        assert round(float(protein_stat_res[1].pvalue), 4) == 1.0
        assert protein_corrected == ["nan", 1.0, "nan", "nan"]

    def test_calculate_compositional_bias_per_site_ascii_path_avoids_unicode_isin(self, mocker, args):
        service = CompositionalBiasPerSite(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AnXT"), id="t1"),
                SeqRecord(Seq("ACxT"), id="t2"),
                SeqRecord(Seq("AG-T"), id="t3"),
            ]
        )
        mocker.patch(
            "phykit.services.alignment.compositional_bias_per_site.np.isin",
            side_effect=AssertionError("ASCII path should use the lookup table"),
        )

        stat_res, corrected = service.calculate_compositional_bias_per_site(
            alignment,
            is_protein=False,
        )

        assert [float(result.statistic) for result in stat_res] == [0.0, 0.0, 0.0, 0.0]
        assert corrected == ["nan", 1.0, "nan", "nan"]

    def test_calculate_compositional_bias_per_site_protein_no_gap_skips_valid_mask(
        self, mocker, args
    ):
        service = CompositionalBiasPerSite(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACDEFGHIK"), id="t1"),
                SeqRecord(Seq("LMNPQRSTV"), id="t2"),
                SeqRecord(Seq("WYACDEFGH"), id="t3"),
            ]
        )
        spy = mocker.spy(cbps_module, "_column_count_stats_from_ascii_codes")

        stat_res, corrected = service.calculate_compositional_bias_per_site(
            alignment,
            is_protein=True,
        )

        assert spy.call_args.args[1] is None
        assert len(stat_res) == 9
        assert len(corrected) == 9

    def test_column_sum_squares_matches_explicit_reduction(self):
        counts = np.array(
            [
                [2, 0, 3, 1],
                [1, 4, 0, 2],
                [0, 1, 5, 3],
            ],
            dtype=np.int64,
        )

        np.testing.assert_allclose(
            cbps_module._column_sum_squares(counts),
            np.sum(counts * counts, axis=0),
        )

    def test_column_totals_use_array_reduction_for_common_counts(self, monkeypatch):
        counts = np.arange(4 * 12000, dtype=np.float64).reshape(4, 12000)
        expected = counts.sum(axis=0)

        def fail_sum(*_args, **_kwargs):
            raise AssertionError("common compositional-bias totals should use ndarray.sum")

        monkeypatch.setattr(cbps_module.np, "sum", fail_sum)

        np.testing.assert_allclose(cbps_module._column_totals(counts), expected)

    def test_column_totals_preserve_long_matrix_sum_path(self, monkeypatch):
        counts = np.ones((4, 20001), dtype=np.float64)
        original_sum = cbps_module.np.sum
        calls = []

        def sum_spy(values, *args, **kwargs):
            calls.append((values.shape, kwargs.get("axis")))
            return original_sum(values, *args, **kwargs)

        monkeypatch.setattr(cbps_module.np, "sum", sum_spy)

        np.testing.assert_allclose(
            cbps_module._column_totals(counts),
            np.full(20001, 4.0),
        )
        assert calls == [((4, 20001), 0)]

    def test_column_count_stats_from_ascii_codes_matches_repeated_counts(self, mocker):
        alignment_array = np.array(
            [
                list(b"ACDEFX"),
                list(b"ACGHIX"),
                list(b"LMNPQ-"),
            ],
            dtype=np.uint8,
        )
        gap_lookup = cbps_module._get_gap_lookup(is_protein=True)
        valid_mask = ~gap_lookup[alignment_array]
        valid_symbols = np.unique(alignment_array[valid_mask])
        mocker.patch.object(
            cbps_module.np,
            "tile",
            side_effect=AssertionError("column count helper should avoid tiled offsets"),
        )

        category_counts, totals, sum_squares = (
            cbps_module._column_count_stats_from_ascii_codes(
                alignment_array,
                valid_mask,
                valid_symbols,
                block_size=2,
            )
        )
        counts = np.array(
            [
                np.sum(alignment_array == symbol, axis=0)
                for symbol in valid_symbols
            ],
            dtype=np.float64,
        )

        np.testing.assert_array_equal(
            category_counts,
            np.count_nonzero(counts, axis=0),
        )
        np.testing.assert_allclose(totals, np.sum(counts, axis=0))
        np.testing.assert_allclose(sum_squares, np.sum(counts * counts, axis=0))

    def test_column_count_stats_from_ascii_codes_counts_all_entries_without_mask(self):
        alignment_array = np.array(
            [
                list(b"ACDEF"),
                list(b"ACGHI"),
                list(b"LMNPQ"),
            ],
            dtype=np.uint8,
        )
        valid_symbols = np.unique(alignment_array)

        category_counts, totals, sum_squares = (
            cbps_module._column_count_stats_from_ascii_codes(
                alignment_array,
                None,
                valid_symbols,
                block_size=2,
            )
        )
        counts = np.array(
            [
                np.sum(alignment_array == symbol, axis=0)
                for symbol in valid_symbols
            ],
            dtype=np.float64,
        )

        np.testing.assert_array_equal(
            category_counts,
            np.count_nonzero(counts, axis=0),
        )
        np.testing.assert_allclose(totals, np.sum(counts, axis=0))
        np.testing.assert_allclose(sum_squares, np.sum(counts * counts, axis=0))

    def test_column_count_stats_from_ascii_codes_default_blocks(
        self, monkeypatch
    ):
        width = cbps_module._ASCII_COUNT_BLOCK_SIZE + 1
        alignment_array = np.tile(
            np.frombuffer(b"ACDEFGHIKLMNPQRSTVWY", dtype=np.uint8),
            (3, width // 20 + 1),
        )[:, :width]
        valid_symbols = np.unique(alignment_array)
        original_bincount = cbps_module.np.bincount
        observed_widths = []

        def bincount_spy(values, *args, **kwargs):
            observed_widths.append(kwargs["minlength"] // 256)
            return original_bincount(values, *args, **kwargs)

        monkeypatch.setattr(cbps_module.np, "bincount", bincount_spy)

        cbps_module._column_count_stats_from_ascii_codes(
            alignment_array,
            None,
            valid_symbols,
        )

        assert observed_widths == [cbps_module._ASCII_COUNT_BLOCK_SIZE, 1]

    def test_calculate_compositional_bias_per_site_protein_uses_chunked_counts(
        self, mocker, args
    ):
        service = CompositionalBiasPerSite(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACDEFGHIK"), id="t1"),
                SeqRecord(Seq("LMNPQRSTV"), id="t2"),
                SeqRecord(Seq("WYACDEFGH"), id="t3"),
            ]
        )
        spy = mocker.spy(cbps_module, "_column_count_stats_from_ascii_codes")

        stat_res, corrected = service.calculate_compositional_bias_per_site(
            alignment,
            is_protein=True,
        )

        assert spy.call_count == 1
        assert len(stat_res) == 9
        assert len(corrected) == 9

    def test_plot_compositional_bias_manhattan_creates_file(self, tmp_path):
        pytest.importorskip("matplotlib")
        output = tmp_path / "cbps_plot.png"
        service = CompositionalBiasPerSite(
            Namespace(alignment="x.fa", plot=True, plot_output=str(output))
        )
        rows = [
            {"site": 1, "chi_square": 1.0, "p_value_corrected": 0.01, "p_value": 0.01},
            {"site": 2, "chi_square": 0.5, "p_value_corrected": 0.5, "p_value": 0.5},
            {"site": 3, "chi_square": 0.2, "p_value_corrected": None, "p_value": None},
        ]
        service._plot_compositional_bias_manhattan(rows)
        assert output.exists()

    def test_plot_compositional_bias_skips_redundant_tight_layout(
        self, tmp_path, monkeypatch
    ):
        pytest.importorskip("matplotlib")
        from matplotlib.figure import Figure

        output = tmp_path / "cbps_plot.png"
        service = CompositionalBiasPerSite(
            Namespace(alignment="x.fa", plot=True, plot_output=str(output))
        )
        rows = [
            {"site": 1, "chi_square": 1.0, "p_value_corrected": 0.01, "p_value": 0.01},
            {"site": 2, "chi_square": 0.5, "p_value_corrected": 0.5, "p_value": 0.5},
        ]

        def fail_tight_layout(self, *args, **kwargs):
            raise AssertionError("bbox_inches='tight' handles saved bounds")

        monkeypatch.setattr(Figure, "tight_layout", fail_tight_layout)
        service._plot_compositional_bias_manhattan(rows)

        assert output.exists()

    def test_plot_compositional_bias_manhattan_importerror(self, monkeypatch, capsys):
        original_import = builtins.__import__

        def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name.startswith("matplotlib"):
                raise ImportError("no matplotlib")
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", fake_import)
        service = CompositionalBiasPerSite(Namespace(alignment="x.fa"))
        with pytest.raises(SystemExit) as exc:
            service._plot_compositional_bias_manhattan([{"site": 1, "chi_square": 0.0, "p_value_corrected": 0.5, "p_value": 0.5}])
        assert exc.value.code == 2
        out, _ = capsys.readouterr()
        assert "matplotlib is required for --plot in compositional_bias_per_site" in out

    def test_run_text_output_with_plot_message(self, mocker, capsys):
        service = CompositionalBiasPerSite(
            Namespace(alignment="x.fa", json=False, plot=True, plot_output="cbps.png")
        )
        mocker.patch.object(
            CompositionalBiasPerSite, "get_alignment_and_format", return_value=(_alignment(), "fasta", False)
        )
        mocker.patch.object(
            CompositionalBiasPerSite,
            "calculate_compositional_bias_per_site",
            return_value=(
                [type("R", (), {"statistic": 2.0, "pvalue": 0.0456})()],
                [0.1234],
            ),
        )
        mocked_plot = mocker.patch.object(
            CompositionalBiasPerSite,
            "_plot_compositional_bias_corrected_pvalues",
        )
        mocker.patch.object(
            service,
            "_build_rows",
            side_effect=AssertionError("non-json plot output should not build row dictionaries"),
        )

        service.run()

        out, _ = capsys.readouterr()
        assert "1\t2.0\t0.1234\t0.0456\n" in out
        assert "Saved compositional bias plot: cbps.png" in out
        mocked_plot.assert_called_once_with([0.1234])
