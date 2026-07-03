from argparse import Namespace
import builtins
import subprocess
import sys

import numpy as np
import pytest

from phykit.services.alignment.plot_alignment_qc import PlotAlignmentQC
import phykit.services.alignment.plot_alignment_qc as plot_alignment_qc_module


def test_module_import_defers_heavy_helpers():
    code = """
import sys
import phykit.services.alignment.plot_alignment_qc as module
assert hasattr(module.np, "__getattr__")
assert callable(module.print_json)
assert callable(module.AlignmentOutlierTaxa)
assert "typing" not in sys.modules
assert "numpy" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
assert "phykit.services.alignment.alignment_outlier_taxa" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


class TestPlotAlignmentQC:
    def test_robust_feature_center_and_sigma_large_equal_skips_median_and_std(
        self, monkeypatch
    ):
        values = np.full(
            plot_alignment_qc_module._ROBUST_SIGMA_EQUAL_CHECK_MIN_SIZE + 1,
            0.75,
        )

        def fail_reduction(*_args, **_kwargs):
            raise AssertionError("large equal features should skip median/std")

        monkeypatch.setattr(plot_alignment_qc_module.np, "median", fail_reduction)
        monkeypatch.setattr(plot_alignment_qc_module.np, "std", fail_reduction)

        median, sigma = plot_alignment_qc_module._robust_feature_center_and_sigma(
            values
        )

        assert median == 0.75
        assert sigma == 0.0

    def test_robust_feature_center_and_sigma_matches_mad_reference(self):
        values = np.array([1.0, 2.0, 2.0, 4.0, 100.0])
        median, sigma = plot_alignment_qc_module._robust_feature_center_and_sigma(
            values
        )

        expected_median = float(np.median(values))
        expected_mad = float(np.median(np.abs(values - expected_median)))

        assert median == expected_median
        assert sigma == 1.4826 * expected_mad

    def test_plot_max_uses_array_reduction_for_plot_sized_arrays(self, monkeypatch):
        values = np.array([0.1, 0.9, 0.4])

        def fail_max(*_args, **_kwargs):
            raise AssertionError("plot-sized maxima should use ndarray.max")

        monkeypatch.setattr(plot_alignment_qc_module.np, "max", fail_max)

        assert plot_alignment_qc_module._plot_max(values) == pytest.approx(0.9)

    def test_plot_max_preserves_large_array_np_max_path(self, monkeypatch):
        values = np.ones(
            plot_alignment_qc_module._PLOT_DIRECT_MAX_LIMIT + 1,
            dtype=np.float64,
        )
        original_max = plot_alignment_qc_module.np.max
        calls = []

        def max_spy(observed, *args, **kwargs):
            calls.append(observed.shape)
            return original_max(observed, *args, **kwargs)

        monkeypatch.setattr(plot_alignment_qc_module.np, "max", max_spy)

        assert plot_alignment_qc_module._plot_max(values) == pytest.approx(1.0)
        assert calls == [(plot_alignment_qc_module._PLOT_DIRECT_MAX_LIMIT + 1,)]

    def test_prepare_plot_arrays_reuses_row_values(self):
        rows = [
            {
                "taxon": "t1",
                "occupancy": 1.0,
                "gap_rate": 0.0,
                "composition_distance": 0.1,
                "long_branch_proxy": 0.2,
                "rcvt": 0.01,
                "entropy_burden": 0.02,
            },
            {
                "taxon": "t2",
                "occupancy": 0.7,
                "gap_rate": 0.3,
                "composition_distance": 0.6,
                "long_branch_proxy": None,
                "rcvt": 0.1,
                "entropy_burden": 0.2,
            },
        ]

        arrays = PlotAlignmentQC._prepare_plot_arrays(
            rows,
            [{"taxon": "t2"}],
            ["gap_rate", "long_branch_proxy"],
        )

        assert arrays["taxa"] == ["t1", "t2"]
        np.testing.assert_allclose(arrays["occupancies"], np.array([1.0, 0.7]))
        np.testing.assert_allclose(arrays["gap_rates"], np.array([0.0, 0.3]))
        np.testing.assert_allclose(
            arrays["composition_distances"],
            np.array([0.1, 0.6]),
        )
        np.testing.assert_allclose(
            arrays["long_branch_proxies"],
            np.array([0.2, np.nan]),
            equal_nan=True,
        )
        np.testing.assert_array_equal(arrays["flagged_mask"], np.array([False, True]))
        assert arrays["flagged_taxa"] == {"t2"}
        assert set(arrays["feature_vals"]) == {"gap_rate", "long_branch_proxy"}
        assert arrays["feature_vals"]["gap_rate"] is arrays["gap_rates"]
        assert arrays["feature_vals"]["long_branch_proxy"] is arrays["long_branch_proxies"]

    def test_prepare_plot_arrays_fills_flagged_mask_in_row_pass(self, monkeypatch):
        rows = [
            {
                "taxon": "t1",
                "occupancy": 1.0,
                "gap_rate": 0.0,
                "composition_distance": 0.1,
                "long_branch_proxy": 0.2,
                "rcvt": 0.01,
                "entropy_burden": 0.02,
            },
            {
                "taxon": "t2",
                "occupancy": 0.7,
                "gap_rate": 0.3,
                "composition_distance": 0.6,
                "long_branch_proxy": None,
                "rcvt": 0.1,
                "entropy_burden": 0.2,
            },
        ]

        def fail_fromiter(*_args, **_kwargs):
            raise AssertionError("flagged mask should be filled in the row pass")

        monkeypatch.setattr(plot_alignment_qc_module.np, "fromiter", fail_fromiter)

        arrays = PlotAlignmentQC._prepare_plot_arrays(
            rows,
            [{"taxon": "t2"}],
            ["gap_rate"],
        )

        np.testing.assert_array_equal(arrays["flagged_mask"], np.array([False, True]))

    def test_prepare_plot_arrays_empty_outliers_prefills_false_mask(self):
        rows = [
            {
                "taxon": "t1",
                "occupancy": 0.9,
                "gap_rate": 0.1,
                "composition_distance": 0.2,
                "long_branch_proxy": 0.3,
                "rcvt": 0.04,
                "entropy_burden": 0.05,
            },
            {
                "taxon": "t2",
                "occupancy": 0.8,
                "gap_rate": 0.2,
                "composition_distance": 0.3,
                "long_branch_proxy": None,
                "rcvt": 0.06,
                "entropy_burden": 0.07,
            },
        ]

        arrays = PlotAlignmentQC._prepare_plot_arrays(
            rows,
            [],
            ["occupancy", "entropy_burden"],
        )

        assert arrays["taxa"] == ["t1", "t2"]
        assert arrays["flagged_taxa"] == set()
        np.testing.assert_array_equal(arrays["flagged_mask"], np.array([False, False]))
        np.testing.assert_allclose(arrays["occupancies"], np.array([0.9, 0.8]))
        np.testing.assert_allclose(
            arrays["long_branch_proxies"],
            np.array([0.3, np.nan]),
            equal_nan=True,
        )
        assert arrays["feature_vals"]["occupancy"] is arrays["occupancies"]
        np.testing.assert_allclose(
            arrays["feature_vals"]["entropy_burden"],
            np.array([0.05, 0.07]),
        )

    def test_prepare_plot_arrays_resolves_nan_once_per_call(self, monkeypatch):
        class CountingNumpy:
            def __init__(self, module):
                self._module = module
                self.nan_accesses = 0

            def __getattr__(self, name):
                if name == "nan":
                    self.nan_accesses += 1
                return getattr(self._module, name)

        rows = [
            {
                "taxon": f"t{idx}",
                "occupancy": 0.9,
                "gap_rate": 0.1,
                "composition_distance": 0.2,
                "long_branch_proxy": None,
                "rcvt": 0.04,
                "entropy_burden": 0.05,
            }
            for idx in range(5)
        ]
        counting_np = CountingNumpy(np)
        monkeypatch.setattr(plot_alignment_qc_module, "np", counting_np)

        arrays = PlotAlignmentQC._prepare_plot_arrays(
            rows,
            [{"taxon": "t1"}],
            ["long_branch_proxy"],
        )

        assert counting_np.nan_accesses == 1
        np.testing.assert_allclose(
            arrays["long_branch_proxies"],
            np.full(5, np.nan),
            equal_nan=True,
        )

    def test_composition_distance_panel_batches_scatter_calls(self):
        class FakeAxes:
            def __init__(self):
                self.scatter_calls = []
                self.text_calls = []

            def scatter(self, x, y, **kwargs):
                self.scatter_calls.append((list(x), list(y), kwargs))

            def text(self, x, y, label, **kwargs):
                self.text_calls.append((x, y, label, kwargs))

            def axvline(self, *args, **kwargs):
                pass

            def axvspan(self, *args, **kwargs):
                pass

            def axhline(self, *args, **kwargs):
                pass

            def axhspan(self, *args, **kwargs):
                pass

            def set_title(self, *args, **kwargs):
                pass

            def set_xlabel(self, *args, **kwargs):
                pass

            def set_ylabel(self, *args, **kwargs):
                pass

            def legend(self, *args, **kwargs):
                pass

        rows = [
            {
                "taxon": "t1",
                "composition_distance": 0.1,
                "long_branch_proxy": 0.2,
            },
            {
                "taxon": "t2",
                "composition_distance": 0.6,
                "long_branch_proxy": None,
            },
            {
                "taxon": "t3",
                "composition_distance": 0.5,
                "long_branch_proxy": 0.4,
            },
        ]
        thresholds = {
            "composition_distance": 0.4,
            "long_branch_proxy": 0.3,
        }
        ax = FakeAxes()

        PlotAlignmentQC._plot_composition_distance_panel(
            ax,
            rows,
            {"t2"},
            thresholds,
            normal_color="#636363",
            flagged_color="#ca0020",
            legend_handles=[],
        )

        assert len(ax.scatter_calls) == 2
        assert ax.scatter_calls[0][0] == [0.1, 0.5]
        assert ax.scatter_calls[1][0] == [0.6]
        assert [call[2] for call in ax.text_calls] == ["t2"]

    def test_flag_colors_uses_ordered_mask(self):
        colors = PlotAlignmentQC._flag_colors(
            np.array([False, True, False, True]),
            "#636363",
            "#ca0020",
        )

        assert colors.tolist() == ["#636363", "#ca0020", "#636363", "#ca0020"]

    def test_init(self):
        args = Namespace(
            alignment="x.fa",
            output="out.png",
            width=14.0,
            height=10.0,
            dpi=300,
            gap_z=3.0,
            composition_z=3.0,
            distance_z=3.0,
            rcvt_z=3.0,
            occupancy_z=3.0,
            entropy_z=3.0,
            json=False,
        )
        svc = PlotAlignmentQC(args)
        assert svc.alignment_file_path == "x.fa"
        assert svc.output == "out.png"
        assert svc.width == 14.0
        assert svc.height == 10.0
        assert svc.dpi == 300
        assert svc.gap_z == 3.0
        assert svc.composition_z == 3.0
        assert svc.distance_z == 3.0
        assert svc.rcvt_z == 3.0
        assert svc.occupancy_z == 3.0
        assert svc.entropy_z == 3.0

    def test_process_args_defaults_json_false(self):
        args = Namespace(
            alignment="x.fa",
            output="out.png",
            width=10.0,
            height=8.0,
            dpi=200,
            gap_z=2.0,
            composition_z=2.0,
            distance_z=2.0,
            rcvt_z=2.0,
            occupancy_z=2.0,
            entropy_z=2.0,
        )
        parsed = PlotAlignmentQC(args).process_args(args)
        assert parsed["json_output"] is False

    def test_get_outlier_result_forwards_thresholds(self, monkeypatch):
        captured = {}

        class FakeOutlierService:
            def __init__(self, ns):
                captured["ns"] = ns

            def calculate_outliers(self, alignment, is_protein):
                captured["alignment"] = alignment
                captured["is_protein"] = is_protein
                return {"rows": [], "outliers": [], "thresholds": {}, "features": []}

        monkeypatch.setattr(plot_alignment_qc_module, "AlignmentOutlierTaxa", FakeOutlierService)
        svc = PlotAlignmentQC(
            Namespace(
                alignment="my.fa",
                output="out.png",
                width=14.0,
                height=10.0,
                dpi=300,
                gap_z=1.1,
                composition_z=1.2,
                distance_z=1.3,
                rcvt_z=1.4,
                occupancy_z=1.5,
                entropy_z=1.6,
                json=False,
            )
        )

        result = svc._get_outlier_result("aln", True)

        assert result["rows"] == []
        assert captured["alignment"] == "aln"
        assert captured["is_protein"] is True
        assert captured["ns"].alignment == "my.fa"
        assert captured["ns"].gap_z == 1.1
        assert captured["ns"].composition_z == 1.2
        assert captured["ns"].distance_z == 1.3
        assert captured["ns"].rcvt_z == 1.4
        assert captured["ns"].occupancy_z == 1.5
        assert captured["ns"].entropy_z == 1.6

    def test_run_importerror_exits(self, monkeypatch, capsys):
        original_import = builtins.__import__

        def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name.startswith("matplotlib"):
                raise ImportError("no matplotlib")
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", fake_import)
        svc = PlotAlignmentQC(
            Namespace(
                alignment="x.fa",
                output="out.png",
                width=10.0,
                height=8.0,
                dpi=100,
                gap_z=2.0,
                composition_z=2.0,
                distance_z=2.0,
                rcvt_z=2.0,
                occupancy_z=2.0,
                entropy_z=2.0,
                json=False,
            )
        )

        with pytest.raises(SystemExit) as exc:
            svc.run()
        assert exc.value.code == 2
        out, _ = capsys.readouterr()
        assert "matplotlib is required for plot_alignment_qc" in out

    def test_run_outputs_json(self, tmp_path, monkeypatch):
        pytest.importorskip("matplotlib")
        output_path = tmp_path / "qc.png"
        captured = {}

        class FakeAlignment:
            def __len__(self):
                return 3

            def get_alignment_length(self):
                return 12

        svc = PlotAlignmentQC(
            Namespace(
                alignment="x.fa",
                output=str(output_path),
                width=8.0,
                height=6.0,
                dpi=100,
                gap_z=3.0,
                composition_z=3.0,
                distance_z=3.0,
                rcvt_z=3.0,
                occupancy_z=3.0,
                entropy_z=3.0,
                json=True,
            )
        )
        monkeypatch.setattr(
            svc,
            "get_alignment_and_format",
            lambda: (FakeAlignment(), None, True),
        )
        monkeypatch.setattr(
            svc,
            "_get_outlier_result",
            lambda _alignment, _is_protein: {
                "rows": [
                    {
                        "taxon": "t1",
                        "occupancy": 1.0,
                        "gap_rate": 0.0,
                        "composition_distance": 0.1,
                        "long_branch_proxy": 0.2,
                        "rcvt": 0.01,
                        "entropy_burden": 0.02,
                    },
                    {
                        "taxon": "t2",
                        "occupancy": 0.7,
                        "gap_rate": 0.3,
                        "composition_distance": 0.6,
                        "long_branch_proxy": None,
                        "rcvt": 0.1,
                        "entropy_burden": 0.2,
                    },
                    {
                        "taxon": "t3",
                        "occupancy": 0.8,
                        "gap_rate": 0.2,
                        "composition_distance": 0.5,
                        "long_branch_proxy": 0.4,
                        "rcvt": 0.05,
                        "entropy_burden": 0.1,
                    },
                ],
                "outliers": [{"taxon": "t2"}],
                "thresholds": {
                    "occupancy": 0.75,
                    "gap_rate": 0.25,
                    "composition_distance": 0.4,
                    "long_branch_proxy": 0.3,
                },
                "features": [
                    "gap_rate",
                    "occupancy",
                    "composition_distance",
                    "long_branch_proxy",
                    "rcvt",
                    "entropy_burden",
                ],
            },
        )
        monkeypatch.setattr(plot_alignment_qc_module, "print_json", lambda payload: captured.setdefault("payload", payload))

        svc.run()

        assert output_path.exists()
        assert captured["payload"]["output_file"] == str(output_path)
        assert captured["payload"]["n_taxa"] == 3
        assert captured["payload"]["alignment_length"] == 12
        assert captured["payload"]["outlier_count"] == 1

    def test_run_outputs_terminal_message(self, tmp_path, monkeypatch, capsys):
        pytest.importorskip("matplotlib")
        output_path = tmp_path / "qc.png"

        class FakeAlignment:
            def __len__(self):
                return 1

            def get_alignment_length(self):
                return 4

        svc = PlotAlignmentQC(
            Namespace(
                alignment="x.fa",
                output=str(output_path),
                width=7.0,
                height=5.0,
                dpi=80,
                gap_z=2.0,
                composition_z=2.0,
                distance_z=2.0,
                rcvt_z=2.0,
                occupancy_z=2.0,
                entropy_z=2.0,
                json=False,
            )
        )
        monkeypatch.setattr(svc, "get_alignment_and_format", lambda: (FakeAlignment(), None, True))
        monkeypatch.setattr(
            svc,
            "_get_outlier_result",
            lambda _alignment, _is_protein: {
                "rows": [
                    {
                        "taxon": "t1",
                        "occupancy": 1.0,
                        "gap_rate": 0.0,
                        "composition_distance": 0.1,
                        "long_branch_proxy": 0.2,
                        "rcvt": 0.01,
                        "entropy_burden": 0.02,
                    }
                ],
                "outliers": [],
                "thresholds": {
                    "occupancy": None,
                    "gap_rate": None,
                    "composition_distance": None,
                    "long_branch_proxy": None,
                },
                "features": [
                    "gap_rate",
                    "occupancy",
                    "composition_distance",
                    "long_branch_proxy",
                    "rcvt",
                    "entropy_burden",
                ],
            },
        )

        svc.run()

        out, _ = capsys.readouterr()
        assert output_path.exists()
        assert "Saved alignment QC plot:" in out
